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
        parser = argparse.ArgumentParser(description='Run OrthoMCL on a group of proteomes.',
                                         epilog="""
Stages: 
1 - Preprocessing ('orthomclAdjustFasta' and 'orthomclFilterFasta'), 
2 - BLAST search, 
3 - Database processing, 
4 - Post processing (MCL)
""")
        parser.add_argument("-i", "--input", metavar='<input_table.tsv>', required=True,
                            help="A table containing paths of proteome FASTA in the first column "
                                 "and species abbreviations to use in the second column. "
                                 "The each abbreviation must have length of 4 characters, with capitalized "
                                 "only the first letter, e.g. Ecol (Escherichia coli) or Klsp (Klebsiella sp.)")
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
        self.input_table = self._namespace.input
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


class OrthoMCLHandler:
    def __init__(self):
        self.output_dir = validator.output_dir
        self.sampledata_dict = self._parse_abbreviations_table(validator.input_table)
        self.orthomcl_cfg_file = validator.config_file
        self.sql_keeper = self._parse_db_properties()

    @staticmethod
    def _parse_abbreviations_table(table_file: str):
        logging.info("Parse abbreviations table")
        out = []
        lines = Utils.load_list(table_file)
        for line in lines:
            if len(line.strip()) == 0:
                continue
            pfasta, abbr = [i.strip() for i in line.split()]
            if len(abbr) > 4:  # OrthoMCL requirement
                logging.warning("The abbreviation is too long, it is recommend to truncate it "
                                "to 4 characters: '{}'".format(abbr))
            if len(abbr) < 3:
                Utils.log_and_raise(
                    "The abbreviation is too small, must contain at least 3 characters: '{}'".format(abbr))
            # abbr = abbr[:4]
            abbr = abbr.capitalize()
            out.append({"pfasta": pfasta, "abbr": abbr})
        return out

    @staticmethod
    def _parse_orthomcl_cfg(cfg_file: str):
        import configparser
        logging.info("Parse OrthoMCL configuration file")
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

    @staticmethod
    def _orthomcl_adjust_fasta(d: dict):
        Utils.run_and_log("orthomclAdjustFasta {} {} 1".format(d["abbr"], d["pfasta"]))

    def run_orthomcl_adjust_fasta(self):
        out_dir = os.path.join(self.output_dir, "compliantFasta")
        os.makedirs(out_dir, exist_ok=True)
        os.chdir(out_dir)
        logging.info("Adjust FASTA")
        _ = Utils.single_core_queue(self._orthomcl_adjust_fasta, self.sampledata_dict)
        os.chdir(validator.output_dir)
        logging.info("Filter FASTA")
        Utils.run_and_log("orthomclFilterFasta compliantFasta 10 20")

    def run_diamond(self):
        os.chdir(self.output_dir)
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
        os.chdir(self.output_dir)
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

    def run_mcl(self):
        os.chdir(self.output_dir)
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

    def convert_mcl_groups_to_pivot(self, annotation_file: str, sample_name: str):
        os.chdir(self.output_dir)
        mcl_groups_df = pd.DataFrame(self._parse_mcl_groups("mcl_groups.txt"))
        annotation_df = pd.read_table(annotation_file, encoding="utf-8", sep="\t", header=0)
        annotation_df["sample_name"] = sample_name
        pass

    def fix_permissions(self):
        logging.info("Fix permissions for the output dir: {}".format(self.output_dir))
        Utils.run_and_log("chmod -R 777 {}".format(self.output_dir))

    def handle(self):
        functions = (self.run_orthomcl_adjust_fasta, self.run_diamond, self.run_mysql_tasks, self.run_mcl)
        for idx in validator.stages_to_do:
            logging.info("Start the pipeline step {}".format(idx))
            functions[idx]()
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
