#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
import subprocess
import multiprocessing


class ArgValidator:
    abbreviations_table, stages_to_do, threads_number, config_file, output_dir = (None,) * 5

    def __init__(self):
        parser = argparse.ArgumentParser(description='Run OrthoMCL on a group of proteomes.',
                                         epilog="""
Stages: 
1 - Preprocessing ('orthomclAdjustFasta' and 'orthomclFilterFasta'), 
2 - BLAST search, 
3 - Database processing, 
4 - Post processing (MCL)
""")
        parser.add_argument('-a', '--abbr_table', metavar='<input_table.tsv>', required=True,
                            help="A table containing paths of proteome FASTA in the first column "
                                 "and species abbreviations to use (column 2)")
        parser.add_argument('-s', '--start', help='Stage to start the pipeline', type=int, default=1,
                            metavar='<1|2|3|4>', choices=[1, 2, 3, 4])
        parser.add_argument('-f', '--finish', help='Stage to finish the pipeline', type=int, default=4,
                            metavar='<1|2|3|4>', choices=[1, 2, 3, 4])
        parser.add_argument('-t', '--threads', help='Number of threads to run at the BLAST stage',
                            metavar='<int>', type=int, default=multiprocessing.cpu_count())
        parser.add_argument('-c', '--config', help='The OrthoMCL config file', metavar='<file>',
                            default="/opt/my_tools/orthomcl.config")
        parser.add_argument('-o', '--output_dir', help='Output directory', metavar='<dir>', required=True,)
        self._namespace = parser.parse_args()
        self.verify()
        os.chdir(self.output_dir)

    @staticmethod
    def log_and_raise(msg):
        logging.CRITICAL(msg)
        raise ValueError(msg)

    def verify(self):
        self.abbreviations_table = self._namespace.table
        start_point = self._namespace.start
        finish_point = self._namespace.finish
        if start_point > finish_point:
            self.log_and_raise("Start stage ({}) must be before end stage ({})".format(start_point, finish_point))
        # Make 'stages_to_do' zero-based
        self.stages_to_do = range(start_point - 1, finish_point)
        self.threads_number = self._namespace.threads
        total_threads = multiprocessing.cpu_count()
        if self.threads_number > total_threads:
            self.threads_number = total_threads
            logging.WARNING("The given threads number ({}) is too large, using {} threads by default".format(
                self.threads_number, total_threads))
        self.config_file = self._namespace.config
        if not os.path.isfile(self.config_file):
            self.log_and_raise("Not found: '{}'".format(self.config_file))
        self.output_dir = os.path.normpath(self.output_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        else:
            logging.WARNING("The output directory exists: '{}'".format(self.output_dir))


class OrthoMCLHandler:
    @staticmethod
    def single_core_queue(func, queue):
        return [func(i) for i in queue]

    @staticmethod
    def log_and_process(cmd: str, log_file: str = None):
        logging.DEBUG("Processing command '{}'".format(cmd))
        log = subprocess.getoutput(cmd)
        if not log_file:
            logging.debug(log)
        else:
            with open(log_file, mode="w", encoding="utf-8") as f:
                f.write(log)
                f.close()

    @staticmethod
    def _parse_abbreviations_table(table_file):
        out = []
        with open(table_file, mode="r", encoding="utf-8") as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                d = dict()
                d["pfasta"], d["abbr"] = [i.strip() for i in line.split()]
                out.append(d)
            f.close()
        return out

    @staticmethod
    def _orthomcl_adjust_fasta(d):
        OrthoMCLHandler.log_and_process("orthomclAdjustFasta '{}' '{}' 1".format(d["abbr"], d["pfasta"]))

    def run_orthomcl_adjust_fasta(self):
        out_dir = os.path.join(validator.output_dir, "compliantFasta")
        os.makedirs(out_dir, exist_ok=True)
        os.chdir(out_dir)
        logging.INFO("Parse abbreviations table")
        queue = self._parse_abbreviations_table(validator.abbreviations_table)
        logging.INFO("Adjust FASTA")
        _ = self.single_core_queue(self._orthomcl_adjust_fasta, queue)
        os.chdir(validator.output_dir)
        logging.INFO("Filter FASTA")
        self.log_and_process("orthomclFilterFasta compliantFasta 10 20")

    def run_diamond(self):
        os.chdir(validator.output_dir)
        logging.INFO("Make diamond database")
        self.log_and_process("diamond makedb --in goodProteins.fasta -d goodProteins.fasta")
        logging.INFO("Do diamond search")
        self.log_and_process("diamond blastp --threads {} \
        --db goodProteins.fasta --outfmt 6 --out all_v_all.blastp \
        --query goodProteins.fasta --max-target-seqs 100000 --evalue 1e-5 --masking 1".format(validator.threads_number))
        logging.INFO("Process blast results to the MySQL-ready file")
        self.log_and_process("orthomclBlastParser all_v_all.blastp compliantFasta >> similarSequences.txt")

    @staticmethod
    def _parse_orthomcl_cfg(cfg_file):
        import configparser
        with open(cfg_file, mode="r", encoding="utf-8") as f:
            # configparser is sensitive to header
            cfg_buf = "\n".join(["[Main]"] + [j for j in [i.strip() for i in f] if len(j) > 0])
            f.close()
        cp = configparser.ConfigParser()
        cp.read_string(cfg_buf)
        out = {i[0]: i[1] for i in cp.items("Main")}
        return out

    def run_mysql_tasks(self):
        os.chdir(validator.output_dir)
        _CFG = validator.config_file
        cfg_dict = self._parse_orthomcl_cfg(_CFG)
        _HOST = "localhost"
        _DB = cfg_dict.get("dbconnectstring").split(":")[2]
        _USER = cfg_dict.get("dblogin")
        _PASSWORD = cfg_dict.get("dbpassword")
        logging.INFO("Delete the old MySQL database")
        self.log_and_process("mysql --host {} -u {} -p{} -e 'DROP DATABASE IF EXISTS {}'".format(_HOST, _USER,
                                                                                                 _PASSWORD, _DB))
        logging.INFO("Create the new MySQL database")
        self.log_and_process("mysql --host {} -u {} -p{} -e 'CREATE DATABASE {}'".format(_HOST, _USER, _PASSWORD, _DB))
        logging.INFO("Install OrthoMCL schema")
        self.log_and_process("orthomclInstallSchema {} orthomclInstallSchema.log".format(_CFG))
        logging.INFO("Push data into the database")
        self.log_and_process("orthomclLoadBlast {} similarSequences.txt".format(_CFG))
        logging.INFO("Process pairs")
        self.log_and_process("orthomclPairs {} orthomclPairs.log cleanup=yes".format(_CFG))
        logging.INFO("Dump pairs")
        self.log_and_process("orthomclDumpPairsFiles {}".format(_CFG))

    def run_mcl(self):
        os.chdir(validator.output_dir)
        logging.INFO("Run MCL")
        self.log_and_process("mcl mclInput --abc -I 1.5 -o mclOutput")
        logging.info("Make groups file")
        self.log_and_process("orthomclMclToGroups groups 1000 < mclOutput > groups.txt")

    def fix_permissions(self):
        logging.info("Fix permissions for the output dir: {}".format(validator.output_dir))
        self.log_and_process("chmod -R 777 {}".format(validator.output_dir))

    def handle(self):
        functions = (self.run_orthomcl_adjust_fasta, self.run_diamond, self.run_mysql_tasks, self.run_mcl)
        for idx in validator.stages_to_do:
            logging.INFO("Start the pipeline step {}".format(idx))
            functions[idx]()
        self.fix_permissions()


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    validator = ArgValidator()
    handler = OrthoMCLHandler()
    handler.handle()
    logging.INFO("The pipeline processing has been completed")
