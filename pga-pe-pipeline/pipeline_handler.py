#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOT to be launched with Docker

import os
import sys
import logging
import subprocess


class SampleDataLine:
    def __init__(self, sample_name: str, sample_reads: list, taxa: list):
        # e.g "ecoli_sample", ["reads.1.fq", "reads.2.fq"], ["Escherichia", "coli", "O157:H7"]
        self.sample_name = sample_name
        self.sample_reads = sample_reads
        self.taxa = taxa
        _MAX_PREFIX_LENGTH = 4
        if len(self.taxa[1]) >= _MAX_PREFIX_LENGTH - 1:
            self.prefix = self.taxa[0][0].upper() + self.taxa[1][:_MAX_PREFIX_LENGTH - 1].lower()
        else:
            self.prefix = self.taxa[0][:_MAX_PREFIX_LENGTH - len(self.taxa[1])].capitalize() + self.taxa[1].lower()
    def update_reads(self, reads: list):
        return SampleDataLine(sample_name=self.sample_name, sample_reads=reads, taxa=self.taxa.copy())


class Handler:
    TOOLS = ("FastQC", "Trimmomatic", "cutadapt", "SPAdes", "Prokka", "bowtie2", "SnpEff", "srst2", "OrthoMCL")
    DOCKER_RUN_CMD = "docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it"
    def __init__(self, output_dir: str):
        self.output_dir_root = os.path.normpath(output_dir)
        if os.path.exists(self.output_dir_root):
            logging.WARNING("The path exists: '{}'".format(self.output_dir_root))
        os.makedirs(self.output_dir_root, exist_ok=True)
        # Output paths for each step
        self.output_dirs = {i: os.path.normpath(os.path.join(self.output_dir_root, "{}_{}".format(idx, i))) for idx, i in self.TOOLS}
    @staticmethod
    def filename_only(path):
        return os.path.splitext(os.path.basename(os.path.normpath(path)))[0]
    @staticmethod
    def get_page(url, tries: int = 5):
        def _get_page(_url):
            return subprocess.getoutput("curl -fsSL '{}'".format(url))
        out = _get_page(url)
        _try = 0
        while out.startswith("curl") and _try < tries:  # An error report
            out = _get_page(url)
            _try += 1
        return out
    def run_quay_image(self, img_name, img_tag: str = None, repo_name: str = "biocontainers", cmd: str = "echo"):
        import json
        if not img_tag:
            api_response = json.loads(self.get_page("https://quay.io/api/v1/repository/{}/{}".format(repo_name, img_name)))
            img_tag = sorted(set(api_response.get("tags")))[-1]
        img_name_full = "quay.io/{}/{}:{}".format(repo_name, img_name, img_tag)
        print("Using image: '{}'".format(img_name_full))
        docker_cmd = "docker pull {a} && {b} {a}".format(a=img_name_full, b=self.DOCKER_RUN_CMD)
        out_cmd = docker_cmd.strip() + " " + cmd.strip()  # cmd may contain curly braces, so str.format() is not usable
        print(out_cmd)
        return subprocess.getoutput(out_cmd)
    @staticmethod
    def clean_path(path):
        if os.path.exists(path):
            logging.WARNING("The path exists and will be replaced with the new data: '{}'".format(path))
            if path != "/":  # I know this is useless
                _ = subprocess.getoutput("rm -rf {}".format(os.path.normpath(path)))
        os.makedirs(path, exist_ok=True)  # path could be created under root
    # Pipeline steps
    def run_fastqc(self, single_reads_file: str):
        # One per read file (two per paired-end sample)
        stage_dir = os.path.join(self.output_dirs["FastQC"], self.filename_only(single_reads_file))
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        cmd = 'bash -c "fastqc {} && chmod -R 777 {}"'.format(single_reads_file, stage_dir)
        os.chdir(self.output_dir_root)
        log = self.run_quay_image("fastqc", cmd=cmd)
        print(log)
    def run_trimmomatic(self, sampledata: SampleDataLine):
        # One per sample
        stage_dir = os.path.join(self.output_dirs["Trimmomatic"])
        self.clean_path(stage_dir)
    def run_prokka(self, genus: str, species: str, sample_name: str, assembly: str, out_dir: str):
        out_dir = "\'{}\'".format(os.path.normpath(out_dir))
        cmd = """
        bash -c \
            "mkdir -p {o} && \
             cd {o} && \
             prokka --cpu $(nproc) --outdir {o} --force --prefix {n} --locustag \'prokka\' \
                    --genus {g} --species {s} \'{a}\' && \
             chmod -R 777 {o}"
        """.format(a=assembly, n=sample_name, g=genus, s=species, o=out_dir)
        log = self.run_quay_image("prokka", cmd=cmd)
        print(log)
    def run_srst2(self, genus: str, species: str, pe_fasta_1: str, pe_fasta_2: str, out_dir: str):
        genus = genus.strip().capitalize()
        species = species.strip().lower()
        out_dir = "\'{}\'".format(os.path.normpath(out_dir))
        mlst_db_local = "\'{}_{}.fasta\'".format(genus, species)
        mlst_definitions_local = "\'{}{}.txt\'".format(genus.lower()[0], species)
        # Default cmd
        srst2_cmd_processed = """
        srst2 --output test --input_pe *.fastq.gz --mlst_db {db} --mlst_definitions {de} --mlst_delimiter '_'
        """.format(db=mlst_db_local, de=mlst_definitions_local)
        if not all(
                os.path.isfile(i) for i in (os.path.join(out_dir, i) for i in (mlst_db_local, mlst_definitions_local))):
            getmlst_cmd = """
            bash -c \
                "mkdir -p {o} && \
                 cd {o} && \
                 getmlst.py --species \'{g} {s}\'"
            """.format(g=genus, s=species, o=out_dir)
            getmlst_log = self.run_quay_image("srst2", cmd=getmlst_cmd)
            srst2_cmd_raw = getmlst_log.split("Suggested srst2 command for use with this MLST database:")[-1].strip()
            if not srst2_cmd_raw.startswith("srst2"):
                print(getmlst_log)
                raise ValueError("`getmlst.py` did not finished correctly!")
            srst2_cmd_processed = srst2_cmd_raw.replace("*.fastq.gz", " ".join([pe_fasta_1, pe_fasta_2])).replace(
                "--output test", "--output {}_MLST".format(genus))
        srst2_cmd_processed = """
        bash -c \
            "cd {o} && \
             {c} --forward _R1_001 --reverse _R2_001 --log && \
             chmod -R a+rw {o}"
        """.format(c=srst2_cmd_processed, o=out_dir)
        srst2_log = self.run_quay_image("srst2", cmd=srst2_cmd_processed)
        print(srst2_log)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.INFO("The pipeline processing has been started")
    logging.INFO("The pipeline processing has been completed")
