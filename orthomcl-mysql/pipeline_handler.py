#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import argparse
import logging
import subprocess
import multiprocessing
import pandas as pd
import pymysql


class ArgValidator:
    def __init__(self):
        parser = argparse.ArgumentParser(description="Run OrthoMCL on a group of proteomes.",
                                         epilog="""
Stages: 
1 - Preprocessing ('orthomclAdjustFasta' and 'orthomclFilterFasta'), 
2 - BLAST search, 
3 - Database processing, 
4 - Post processing (MCL)
""")
        parser.add_argument("-i", "--input", metavar="<input.sampledata>", required=True,
                            help="A tab-delimited table containing two columns of data for each sample per row. "
                                 "The first column must contain unique sample number, strain name or else identifier. "
                                 "The second column is intended for paths of annotated GenBank data files.")
        parser.add_argument('-s', '--start', help='Stage to start the pipeline', type=int, default=1,
                            metavar='<1|2|3|4>', choices=[1, 2, 3, 4])
        parser.add_argument('-f', '--finish', help='Stage to finish the pipeline', type=int, default=4,
                            metavar='<1|2|3|4>', choices=[1, 2, 3, 4])
        parser.add_argument('-t', '--threads', help='Number of threads to run at the BLAST stage',
                            metavar='<int>', type=int, default=multiprocessing.cpu_count())
        parser.add_argument('-c', '--config', help='The OrthoMCL config file', metavar='<file>',
                            default="/opt/my_tools/orthomcl.config")
        parser.add_argument("-o", "--output_dir", help='Output directory', metavar='<dir>', required=True)
        self._namespace = parser.parse_args()
        self.input_sampledata = self._namespace.input
        self.stages_to_do = []
        self.threads_number = self._namespace.threads
        self.config_file = self._namespace.config
        self.output_dir = os.path.normpath(self._namespace.output_dir)
        self.verify()
        os.chdir(self.output_dir)

    def verify(self):
        start_point = self._namespace.start
        finish_point = self._namespace.finish
        if start_point > finish_point:
            Utils.log_and_raise("Start stage ({}) must be before end stage ({})".format(start_point, finish_point))
        # Make 'stages_to_do' zero-based
        self.stages_to_do = range(start_point - 1, finish_point)
        total_threads = multiprocessing.cpu_count()
        if self.threads_number > total_threads:
            self.threads_number = total_threads
            logging.warning("The given threads number ({}) is too large, using {} threads by default".format(
                self.threads_number, total_threads))
        if not os.path.isfile(self.config_file):
            Utils.log_and_raise("Not found: '{}'".format(self.config_file))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        else:
            logging.warning("The output directory exists: '{}'".format(self.output_dir))


class SQLPropertiesKeeper:
    host, db, user, password = ("",) * 4


class MySQLQueryManager:
    def __init__(self, keeper: SQLPropertiesKeeper):
        self.keeper = keeper

    def execute(self, queries: list):
        con = pymysql.connect(host=self.keeper.host, user=self.keeper.user, password=self.keeper.password,
                              db=self.keeper.db, charset="utf8mb4", cursorclass=pymysql.cursors.DictCursor)
        with con.cursor() as cur:
            for query in queries:
                cur.execute(query)
            out = [i for i in cur]
        con.commit()
        con.close()
        return out


class SampleDataLine:
    def __init__(self, name: str, genbank: str):
        self.name = name
        self.genbank = genbank
        self.pfasta = ""
        self.annotation = ""
        self.is_valid = True
        self.validate()

    def validate(self):
        if not any(os.path.isfile(i) for i in (self.genbank, self.pfasta)):
            logging.warning("Not found any of the data files: GenBank: '{}', protein FASTA: '{}'".format(
                self.genbank, self.pfasta))
            self.is_valid = False

    @staticmethod
    def parse(sampledata_row: str):
        name, genbank = Utils.remove_empty_values(sampledata_row.strip().split("\t"))
        if len(name) > 4:  # OrthoMCL requirement
            logging.warning("The abbreviation is too long, it is recommend to truncate it "
                            "to 4 characters: '{}'".format(name))
        if len(name) < 3:
            logging.warning(
                "The abbreviation is too small, it is recommended to contain at least 3 characters: '{}'".format(name))
        return SampleDataLine(name, genbank)


class SampleDataArray:
    lines = []

    def __len__(self):
        return len(self.lines)

    def get_names(self):
        return [i.name for i in self.lines]

    def validate(self):
        self.lines = [i for i in self.lines if i.is_valid]
        if len(self.lines) == 0:
            Utils.log_and_raise("No valid protein FASTA sequence files were found, exit.")
        names = self.get_names()
        for name in set(names):
            if names.count(name) > 1:
                Utils.log_and_raise("Duplicate sample name found, please rename it: '{}'".format(name))
        logging.info("{} input files will be processed".format(len(self)))

    def apply_single_core_function(self, func):
        _ = Utils.single_core_queue(func, queue=self.lines)

    @staticmethod
    def parse(sampledata_file: str):
        logging.info("Parse input sample data: '{}'".format(sampledata_file))
        all_lines = Utils.remove_empty_values(Utils.load_list(sampledata_file))
        array = SampleDataArray()
        array.lines = [SampleDataLine.parse(i) for i in all_lines]
        array.validate()
        return array


class OrthoMCLHandler:
    def __init__(self):
        self.output_dir_root = validator.output_dir
        self.sampledata_array = SampleDataArray.parse(validator.input_sampledata)
        self.orthomcl_cfg_file = validator.config_file
        self.sql_keeper = self._parse_db_properties()

    @staticmethod
    def _parse_orthomcl_cfg(cfg_file: str):
        import configparser
        logging.info("Parse OrthoMCL configuration file: '{}'".format(cfg_file))
        # configparser is sensitive to header
        cfg_buf = "\n".join(["[Main]"] + Utils.load_list(cfg_file))
        cp = configparser.ConfigParser()
        cp.read_string(cfg_buf)
        out = {i[0]: i[1] for i in cp.items("Main")}
        return out

    def _parse_db_properties(self):
        db_config_dict = self._parse_orthomcl_cfg(self.orthomcl_cfg_file)
        keeper = SQLPropertiesKeeper()
        # dbConnectString=dbi:mysql:orthomcl:localhost:3306
        db_connect_list = db_config_dict.get("dbconnectstring").split(":")
        keeper.db = db_connect_list[2]
        keeper.host = db_connect_list[3]
        keeper.port = db_connect_list[4]
        keeper.user = db_config_dict.get("dblogin")
        keeper.password = db_config_dict.get("dbpassword")
        return keeper

    def _extract_pfasta_from_gbk(self, sampledata: SampleDataLine):
        _TOOL = "extract_pfasta_from_gbk"
        tool_dir = os.path.join(self.output_dir_root, _TOOL)
        os.makedirs(tool_dir, exist_ok=True)
        sampledata.pfasta = os.path.join(tool_dir, "{}.faa".format(sampledata.name))
        # The following directive is actually taken from the corresponding script
        sampledata.annotation = "{}_annotation.tsv".format(".".join(sampledata.pfasta.split(".")[:-1]))
        sampledata.annotation = os.path.join(self.output_dir_root, sampledata.annotation)
        logging.debug("Extract amino acid FASTA for the GenBank data file: '{}'".format(sampledata.genbank))
        Utils.run_and_log("python3 /opt/my_tools/{}.py -i {} -s {} -o {}".format(
            _TOOL, sampledata.genbank, sampledata.name, sampledata.pfasta))

    def extract_pfasta_from_gbk_wrapper(self):
        logging.info("Create protein sequence and annotation files")
        self.sampledata_array.apply_single_core_function(self._extract_pfasta_from_gbk)

    @staticmethod
    def _orthomcl_adjust_fasta(line: SampleDataLine):
        Utils.run_and_log("orthomclAdjustFasta {} {} 1".format(line.name, line.pfasta))

    def run_orthomcl_adjust_fasta(self):
        tool_dir = os.path.join(self.output_dir_root, "compliantFasta")
        os.makedirs(tool_dir, exist_ok=True)
        os.chdir(tool_dir)
        logging.info("Adjust FASTA")
        self.sampledata_array.apply_single_core_function(self._orthomcl_adjust_fasta)
        os.chdir(self.output_dir_root)
        logging.info("Filter FASTA")
        Utils.run_and_log("orthomclFilterFasta compliantFasta 10 20")

    def run_diamond(self):
        os.chdir(self.output_dir_root)
        logging.info("Make diamond database")
        Utils.run_and_log("diamond makedb --in goodProteins.fasta -d goodProteins.fasta")
        logging.info("Do diamond search")
        Utils.run_and_log("""
        diamond blastp --threads {} \
            --db goodProteins.fasta --outfmt 6 --out all_v_all.blastp \
            --query goodProteins.fasta --max-target-seqs 100000 --evalue 1e-5 --masking 1
        """.format(validator.threads_number))
        logging.info("Process blast results to the MySQL-ready file")
        Utils.run_and_log("orthomclBlastParser all_v_all.blastp compliantFasta >> similarSequences.txt")

    def run_mysql_tasks(self):
        os.chdir(self.output_dir_root)
        logging.info("Check the MySQL port")
        Utils.run_and_log("""mysql -e 'SHOW GLOBAL VARIABLES LIKE "PORT"' | grep -Po '[0-9]{2,}'""")
        logging.info("Delete the old MySQL database")
        Utils.run_and_log("mysql --host {} -u {} -p{} -e 'DROP DATABASE IF EXISTS {}'".format(
            self.sql_keeper.host, self.sql_keeper.user, self.sql_keeper.password, self.sql_keeper.db))
        logging.info("Create the new MySQL database")
        Utils.run_and_log("mysql --host {} -u {} -p{} -e 'CREATE DATABASE {}'".format(
            self.sql_keeper.host, self.sql_keeper.user, self.sql_keeper.password, self.sql_keeper.db))
        logging.info("Install OrthoMCL schema")
        Utils.run_and_log("orthomclInstallSchema {} orthomclInstallSchema.log".format(self.orthomcl_cfg_file))
        logging.info("Push data into the database")
        Utils.run_and_log("orthomclLoadBlast {} similarSequences.txt".format(self.orthomcl_cfg_file))
        logging.info("Process pairs")
        Utils.run_and_log("orthomclPairs {} orthomclPairs.log cleanup=yes".format(self.orthomcl_cfg_file))
        logging.info("Dump pairs")
        Utils.run_and_log("orthomclDumpPairsFiles {}".format(self.orthomcl_cfg_file))
        os.chdir(self.output_dir_root)

    def run_mcl(self):
        os.chdir(self.output_dir_root)
        logging.info("Run MCL")
        Utils.run_and_log("mcl mclInput --abc -I 1.5 -o mcl_output.txt")
        logging.info("Convert MCL output file to group IDs file")
        zero_filler = len(str(len(Utils.load_list("mcl_output.txt"))))
        Utils.run_and_log("orthomclMclToGroups MCL_ID_ 1{} < mcl_output.txt > mcl_groups.txt".format("0" * zero_filler))

    @staticmethod
    def _parse_mcl_groups(mcl_groups_file):
        lines = Utils.load_list(mcl_groups_file)
        out = []
        for line in lines:
            line_list = Utils.remove_empty_values(line.split(" "))
            if len(line_list) < 2:
                continue
            mcl_id = line_list[0].replace(":", "")
            pfasta_headers = Utils.remove_empty_values(line_list[1:])
            for pfasta_header in pfasta_headers:
                try:
                    prefix, pfasta_id = pfasta_header.split("|")[:2]
                    out.append({"mcl_id": mcl_id, "sample_name": prefix, "pfasta_id": pfasta_id})
                except IndexError:
                    pass
        return out

    def convert_mcl_groups_to_pivot(self):
        _INDEX_COL_NAMES = ["sample_name", "pfasta_id"]
        os.chdir(self.output_dir_root)
        mcl_groups_ds = pd.DataFrame(self._parse_mcl_groups("mcl_groups.txt")).set_index(_INDEX_COL_NAMES)
        annotation_ds = pd.concat(
            [pd.read_table(i.annotation, encoding="utf-8", sep="\t", header=0) for i in self.sampledata_array.lines],
            axis=0, ignore_index=True, sort=False).set_index(_INDEX_COL_NAMES)
        mcl_annotated_ds = pd.concat([annotation_ds, mcl_groups_ds], axis=1, sort=False)
        mcl_annotated_ds.reset_index().to_csv("mcl_annotated_ds.tsv", encoding="utf-8", sep="\t", index=False,
                                              header=True)

    def fix_permissions(self):
        logging.info("Fix permissions for the output dir: {}".format(self.output_dir_root))
        Utils.run_and_log("chmod -R 777 {}".format(self.output_dir_root))

    def handle(self):
        functions = (self.extract_pfasta_from_gbk_wrapper, self.run_orthomcl_adjust_fasta, self.run_diamond,
                     self.run_mysql_tasks, self.run_mcl, self.convert_mcl_groups_to_pivot)
        # TODO: Apply `for idx in validator.stages_to_do:`
        for idx, func in enumerate(functions):
            logging.info("Start the pipeline step {}".format(idx))
            func()
        self.fix_permissions()


class Utils:
    @staticmethod
    def log_and_raise(msg):
        logging.critical(msg)
        raise ValueError(msg)
    @staticmethod
    def run_and_log(cmd: str, log_file: str = None):
        log_cmd = re.sub("[\r\n ]+", " ", cmd.strip())
        logging.debug("Processing command '{}'".format(log_cmd))
        log = subprocess.getoutput(cmd.strip())
        if len(log.strip()) == 0:
            logging.debug("Done")
        else:
            if not log_file:
                logging.debug(log)
            else:
                Utils.dump_string(string=log, file=log_file)
    @staticmethod
    def single_core_queue(func, queue):
        return [func(i) for i in queue]
    @staticmethod
    def load_string(file: str):
        with open(file=file, mode="r", encoding="utf-8") as f:
            s = f.read()
            f.close()
        return s
    @staticmethod
    def dump_string(string: str, file: str):
        os.makedirs(os.path.dirname(file), exist_ok=True)
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)
            f.close()
    @staticmethod
    def remove_empty_values(input_list: list):
        return [j for j in [i.strip() for i in input_list] if len(j) > 0]
    @staticmethod
    def split_lines(string: str):
        return Utils.remove_empty_values(re.sub("[\r\n]+", "\n", string).split("\n"))
    @staticmethod
    def load_list(file: str):
        return Utils.split_lines(Utils.load_string(file))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
    validator = ArgValidator()
    handler = OrthoMCLHandler()
    handler.handle()
    logging.info("The pipeline processing has been completed")
