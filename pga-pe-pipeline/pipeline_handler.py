#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOT to be launched with Docker

import os
import sys
import logging
import subprocess


class SampleDataLine:
    def __init__(self, sample_name: str, sample_reads: list, taxa: list, genome: str = None, plasmid: str = None):
        # e.g "ecoli_sample", ["reads.1.fq", "reads.2.fq"], ["Escherichia", "coli", "O157:H7"]
        self.name = sample_name.strip()
        self.reads = [i.strip() for i in sample_reads]
        self.extension = self.get_extension(self.reads[0].strip())
        self.taxa = self._parse_taxa(taxa)
        self.prefix = self._set_prefix()
        self.genome = genome
        self.plasmid = plasmid
    @staticmethod
    def get_extension(path):
        import pathlib  # Since Python 3.4
        suf = pathlib.Path(path).suffixes
        if len(suf) > 1:  # e.g. '*.fq.gz'
            return "".join(suf[-2:])
        return "".join(suf)
    @staticmethod
    def _parse_taxa(taxa: list):
        out = {i: "" for i in ("genus", "species", "strain")}
        taxa = [j for j in [i.strip() for i in taxa] if len(j) > 0]
        if len(taxa) == 0:
            raise ValueError("Empty taxon data!")
        out["genus"] = taxa[0].capitalize()
        if len(taxa) > 0:
            if not any(i.isdigit() for i in taxa[1]):
                sp = taxa[1].replace(".", "").lower()
                if sp != "sp":
                    out["species"] = sp
            else:
                out["strain"] = taxa[1]
        if len(taxa) > 1:
            out["strain"] = taxa[2]
        return out
    def _set_prefix(self):
        _MAX_PREFIX_LENGTH = 4
        if len(self.taxa[1]) >= _MAX_PREFIX_LENGTH - 1:
            return self.taxa[0][0].upper() + self.taxa[1][:_MAX_PREFIX_LENGTH - 1].lower()
        else:
            return self.taxa[0][:_MAX_PREFIX_LENGTH - len(self.taxa[1])].capitalize() + self.taxa[1].lower()
    def update_reads(self, reads: list):
        import copy
        dc = copy.deepcopy(self)
        dc.reads = reads
        return dc


class Handler:
    TOOLS = ("FastQC", "Trimmomatic", "cutadapt", "SPAdes", "Prokka", "bowtie2", "SnpEff", "srst2", "OrthoMCL")
    DOCKER_RUN_CMD = "docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it"
    def __init__(self, output_dir: str):
        self.output_dir_root = os.path.normpath(output_dir)
        if os.path.exists(self.output_dir_root):
            logging.WARNING("The path exists: '{}'".format(self.output_dir_root))
        os.makedirs(self.output_dir_root, exist_ok=True)
        # Output paths for each step
        self.output_dirs = {i: os.path.normpath(os.path.join(self.output_dir_root, "{}_{}".format(idx, i))) for idx, i
                            in enumerate(self.TOOLS)}
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
                print(subprocess.getoutput("rm -rf {}".format(os.path.normpath(path))))
        os.makedirs(path, exist_ok=True)  # path could be created under root and could not be deleted easily
    @staticmethod
    def process_reads(sampledata: SampleDataLine, out_dir: str, suffix: str):
        out_dir, suffix = [i.strip() for i in (out_dir, suffix)]
        processed_reads = [
            os.path.join(out_dir, "{}_{}.{}{}".format(sampledata.name, suffix, idx + 1, sampledata.extension)) for
            idx, i in enumerate(sampledata.reads)]
        return sampledata.update_reads(processed_reads)
    # Pipeline steps
    def run_fastqc(self, sampledata: SampleDataLine):
        _TOOL = "fastqc"
        # One per read sample
        for idx, reads_file in enumerate(sampledata.reads):
            stage_dir = os.path.join(self.output_dirs["FastQC"], "{}_{}".format(sampledata.name, idx + 1))
            self.clean_path(stage_dir)
            os.chdir(stage_dir)
            cmd = "bash -c '{} {} && chmod -R 777 {}'".format(_TOOL, reads_file, stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd)
            print(log)
        os.chdir(self.output_dir_root)
        # Reads are unchanged, so there is nothing to return
    def run_trimmomatic(self, sampledata: SampleDataLine):
        # One per sample
        _TOOL = "trimmomatic"
        stage_dir = os.path.join(self.output_dirs["Trimmomatic"], sampledata.name)
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        trimmed_sampledata = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        untrimmed_sampledata = self.process_reads(sampledata, out_dir=stage_dir, suffix="{}_untrinned".format(_TOOL))
        cmd = """
        bash -c \
            '{T} PE -phred33 {r1} {r2} {t1} {u1} {t2} {u2} \
                         ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && \
             chmod -R 777 {o}'
        """.format(T=_TOOL, o=stage_dir, r1=sampledata.reads[0], r2=sampledata.reads[1], t1=trimmed_sampledata.reads[0],
                   t2=trimmed_sampledata.reads[1], u1=untrimmed_sampledata.reads[0], u2=untrimmed_sampledata.reads[1])
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        os.chdir(self.output_dir_root)
        return trimmed_sampledata
    def run_cutadapt(self, sampledata: SampleDataLine):
        # One per sample
        _TOOL = "cutadapt"
        _ADAPTER = "AGATCGGAAGAG"
        stage_dir = os.path.join(self.output_dirs["cutadapt"], sampledata.name)
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        trimmed_sampledata = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        cmd = """
        bash -c 'cutadapt -a {a} -A {a} -m 50 -o {t1} -p {t2} {r1} {r2} && chmod -R 777 {o}'
        """.format(a=_ADAPTER, r1=sampledata.reads[0], r2=sampledata.reads[1], t1=trimmed_sampledata.reads[0],
                   t2=trimmed_sampledata.reads[1], o=stage_dir)
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        os.chdir(self.output_dir_root)
        return trimmed_sampledata
    def run_spades(self, sampledata: SampleDataLine):
        # One per sample
        _TOOL = "spades"
        """
        # Sample launch:
        IMG=quay.io/biocontainers/spades:3.13.1--0 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG bash -c \
            'TOOL=$(find /usr/local/share/ -name spades.py | grep spades | head -n 1) && $TOOL -v'
        """
        stage_dir = os.path.join(self.output_dirs["cutadapt"], sampledata.name)
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        assemblies = {"genome": "", "plasmid": ""}
        for assembly_type in assemblies:
            assembly_dir = os.path.join(stage_dir, assembly_type)
            cmd_append = ""
            if assembly_type == "plasmid":
                cmd_append = " --plasmid"
            cmd = """
            bash -c \
                'TOOL=$(find /usr/local/share/ -name {t}.py | grep {t} | head -n 1) && \
                 $TOOL --careful -o {o} -1 {r1} -2 {r2}{a} && \
                 chmod -R 777 {o}'
            """.format(o=assembly_dir, t=_TOOL, r1=sampledata.reads[0], r2=sampledata.reads[1], a=cmd_append)
            log = self.run_quay_image(_TOOL, cmd=cmd)
            print(log)
            assemblies[assembly_type] = os.path.join(assembly_dir, "contigs.fasta")
        os.chdir(self.output_dir_root)
        sampledata.genome = assemblies["genome"]
        sampledata.plasmid = assemblies["plasmid"]
        return sampledata
    def run_prokka(self, sampledata: SampleDataLine):
        # One per sample
        _TOOL = "prokka"
        """
        # Sample launch:
        IMG=quay.io/biocontainers/prokka:1.12--pl526_0 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG prokka
        """
        stage_dir = os.path.join(self.output_dirs["Prokka"], sampledata.name)
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        taxa_append = ""
        for taxon_name in ("genus", "species", "strain"):
            taxon_value = sampledata.taxa.get(taxon_name)
            if len(taxon_value) > 0:
                taxa_append = "{} --{} {}".format(taxa_append, taxon_name, taxon_value)
        cmd = """
        bash -c \
            'prokka --cpu $(nproc) --outdir {o} --force --prefix {p} --locustag prokka {t} {a} && \
             chmod -R 777 {o}'
        """.format(a=sampledata.genome, t=taxa_append, p=sampledata.prefix, o=stage_dir)
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        os.chdir(self.output_dir_root)
        # TODO Implement return format
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
