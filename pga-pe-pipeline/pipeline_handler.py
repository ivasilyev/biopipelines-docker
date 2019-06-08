#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOT to be launched with Docker
# Get command:
# curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py

import os
import sys
import logging
import subprocess
import multiprocessing


class ArgValidator:
    sampledata_file, output_dir = (None,) * 2
    def __init__(self):
        import argparse
        parser = argparse.ArgumentParser(description="""
Run prokaryotic genome analysis pipeline for group of files with given taxa information""".strip(),
                                         epilog="""
Stages: https://github.com/boulygina/bioinformatics-pipelines/blob/master/Prokaryotes_analysis/prokaryotes.md
            """.strip())
        parser.add_argument(
            "-i", '--input', metavar='<input.sampledata>', required=True,
            help="""
A tab-delimited table file containing information for each sample per row.
Columns:
1. Sample name, e.g. 'eco_01'
2. Paired end reads divided with semicolon, e.g. '/reads/eco_01_R1.fastq.gz;/reads/eco_01_R2.fastq.gz'
3. Taxon information divided with spaces, e.g. 'Escherichia coli O157:H7'""".strip())
        parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Output directory', required=True)
        self._namespace = parser.parse_args()
        self.validate()
    @staticmethod
    def log_and_raise(msg):
        logging.CRITICAL(msg)
        raise ValueError(msg)
    def validate(self):
        self.sampledata_file = self._namespace.input
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        else:
            logging.WARNING("The output directory exists: '{}'".format(self.output_dir))


class SampleDataArray:
    lines = []
    @staticmethod
    def parse(file):
        arr = SampleDataArray()
        with open(file, mode="r", encoding="utf-8") as f:
            arr.lines = [SampleDataLine.parse(j) for j in [i.strip() for i in f] if len(j) > 0]
            f.close()
        return arr


class SampleDataLine:
    prefix, genome, plasmid, annotation_genbank, reference_nfasta, mlst_results = ("", ) * 6
    def __init__(self, sample_name: str, sample_reads: list, taxa: list):
        # e.g "ecoli_sample", ["reads.1.fq", "reads.2.fq"], ["Escherichia", "coli", "O157:H7"]
        self.name = sample_name.strip()
        self.reads = [i.strip() for i in sample_reads]
        self.extension = self.get_extension(self.reads[0].strip())
        self.taxa = self._parse_taxa(taxa)
        self.prefix = self._set_prefix()
    @staticmethod
    def parse(line: str):
        name, reads, taxa = [j for j in [i.strip() for i in line.split("\t")] if len(j) > 0]
        reads = SampleDataLine._parse_reads(reads.split(";"))
        if len(reads) > 2:
            logging.WARNING("More than 2 paired end read files were given: '{}'".format(reads))
        taxa = [j for j in [i.strip() for i in taxa.split(" ")] if len(j) > 0]
        return SampleDataLine(name, reads, taxa)
    @staticmethod
    def get_extension(path):
        import pathlib  # Since Python 3.4
        suf = pathlib.Path(path).suffixes
        if len(suf) > 1:  # e.g. '*.fq.gz'
            return "".join(suf[-2:])
        return "".join(suf)
    @staticmethod
    def _parse_reads(reads: list):
        return [j for j in [i.strip() for i in reads] if len(j) > 0]
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
            self.prefix = self.taxa[0][0].upper() + self.taxa[1][:_MAX_PREFIX_LENGTH - 1].lower()
        else:
            self.prefix = self.taxa[0][:_MAX_PREFIX_LENGTH - len(self.taxa[1])].capitalize() + self.taxa[1].lower()
    def set_reads(self, reads: list):
        reads = self._parse_reads(reads)
        self.reads = reads


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
    def get_page(url, attempts: int = 5):
        def _get_page(_url):
            return subprocess.getoutput("curl -fsSL '{}'".format(url))
        out = _get_page(url)
        _try = 0
        while out.startswith("curl") and _try < attempts:  # An error report
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
        return processed_reads
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
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        untrimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix="{}_untrinned".format(_TOOL))
        cmd = """
        bash -c \
            '{T} PE -phred33 {r1} {r2} {t1} {u1} {t2} {u2} \
                         ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && \
             chmod -R 777 {o}'
        """.format(T=_TOOL, o=stage_dir, r1=sampledata.reads[0], r2=sampledata.reads[1], t1=trimmed_reads[0],
                   t2=trimmed_reads[1], u1=untrimmed_reads[0], u2=untrimmed_reads[1])
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        os.chdir(self.output_dir_root)
        sampledata.set_reads(trimmed_reads)
    def run_cutadapt(self, sampledata: SampleDataLine):
        # One per sample
        _TOOL = "cutadapt"
        _ADAPTER = "AGATCGGAAGAG"
        stage_dir = os.path.join(self.output_dirs["cutadapt"], sampledata.name)
        self.clean_path(stage_dir)
        os.chdir(stage_dir)
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        cmd = """
        bash -c 'cutadapt -a {a} -A {a} -m 50 -o {t1} -p {t2} {r1} {r2} && chmod -R 777 {o}'
        """.format(a=_ADAPTER, r1=sampledata.reads[0], r2=sampledata.reads[1], t1=trimmed_reads[0],
                   t2=trimmed_reads[1], o=stage_dir)
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        os.chdir(self.output_dir_root)
        sampledata.set_reads(trimmed_reads)
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
        taxa_append = ""
        for taxon_name in ("genus", "species", "strain"):
            taxon_value = sampledata.taxa.get(taxon_name)
            if len(taxon_value) > 0:
                taxa_append = "{} --{} {}".format(taxa_append, taxon_name, taxon_value)
        cmd = """
        bash -c \
            'cd {o}
             {t} --cpu {c} --outdir {o} --force --prefix {p} --locustag {p} {a} {g} && chmod -R 777 {o}
             ln -s "{o}/{p}.gbk" "{o}/{p}.gb"
             chmod -R 777 {o}'
        """.format(t=_TOOL, g=sampledata.genome, a=taxa_append, p=sampledata.prefix, o=stage_dir,
                   c=multiprocessing.cpu_count())
        log = self.run_quay_image(_TOOL, cmd=cmd)
        print(log)
        sampledata.annotation_genbank = os.path.join(stage_dir, "{}.gb".format(sampledata.prefix))
    # SNP calling
    def run_bowtie2(self, sampledata: SampleDataLine):
        pass
    def run_samtools(self, sampledata: SampleDataLine):
        pass
    def run_vcftools(self, sampledata: SampleDataLine):
        pass
    # SNP annotation
    def run_snpeff(self, sampledata: SampleDataLine):
        pass
    # MLST typing
    def run_srst2(self, sampledata: SampleDataLine):
        # One per sample, full cleaning is NOT required
        _TOOL = "srst2"
        _GETMLST_ATTEMPTS = 5
        """
        # Sample launch:
        IMG=quay.io/biocontainers/srst2:0.2.0--py27_2 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG getmlst.py -h
        
        docker run --rm --net=host -it $IMG srst2 -h
        """
        tool_dir = self.output_dirs["srst2"]
        stage_dir = os.path.join(tool_dir, sampledata.name)
        self.clean_path(stage_dir)
        genus, species = (sampledata.taxa["genus"], sampledata.taxa["species"])
        mlst_db_local = "{}_{}.fasta".format(genus, species)
        mlst_definitions_local = "{}{}.txt".format(genus.lower()[0], species)
        mlst_db_abs, mlst_definitions_abs = [os.path.join(tool_dir, i) for i in (mlst_db_local, mlst_definitions_local)]
        """
        # Sample output for `getmlst.py --species 'Klebsiella pneumoniae'`
        
        For SRST2, remember to check what separator is being used in this allele database

        Looks like --mlst_delimiter '_'

        >gapA_1  --> -->   ('gapA', '_', '1')

        Suggested srst2 command for use with this MLST database:

        srst2 --output test --input_pe *.fastq.gz --mlst_db Klebsiella_pneumoniae.fasta --mlst_definitions kpneumoniae.txt --mlst_delimiter '_'
        
        # Sample srst2 output:
        
        Klebsiella_MLST__mlst__Klebsiella_pneumoniae__results.txt
        """
        # Default cmd
        srst2_cmd = """
        srst2 --output test --input_pe *.fastq.gz --mlst_db {} --mlst_definitions {} --mlst_delimiter '_'
        """.format(mlst_db_local, mlst_definitions_local)
        if not all(os.path.isfile(i) for i in (mlst_db_abs, mlst_definitions_abs)):
            getmlst_attempt = 0
            while getmlst_attempt < _GETMLST_ATTEMPTS:
                getmlst_attempt += 1
                getmlst_cmd = """
                bash -c \
                    'cd {o} && getmlst.py --species \'{g} {s}\''
                """.format(g=genus, s=species, o=tool_dir)
                getmlst_log = self.run_quay_image("srst2", cmd=getmlst_cmd)
                srst2_cmd = getmlst_log.split("Suggested srst2 command for use with this MLST database:")[-1].strip()
                print(getmlst_log)
                if not srst2_cmd.startswith("srst2"):
                    print("`getmlst.py` has not been finished correctly for attempt {} of {}".format(getmlst_attempt,
                                                                                                     _GETMLST_ATTEMPTS))
        # The input read files must be named by a strict pattern:
        # https://github.com/katholt/srst2#input-read-formats-and-options
        input_reads = ["{}_{}.fastq.gz".format(sampledata.name, idx + 1) for idx, i in enumerate(sampledata.reads)]
        srst2_cmd = srst2_cmd.replace(
            "*.fastq.gz", " ".join(input_reads)).replace(
            "--output test", "--output {}_MLST".format(sampledata.prefix)).replace(
            mlst_db_local, mlst_db_abs).replace(
            mlst_definitions_local, mlst_definitions_abs)
        srst2_cmd_full = """
        bash -c \
            'cd {o} 
             ln -s {r1} {l1}
             ln -s {r2} {l2}
             {c} --log --threads {t}
             chmod -R 777 {o}'
        """.format(c=srst2_cmd.strip(), o=stage_dir, t=multiprocessing.cpu_count(),
                   r1=sampledata.reads[0], r2=sampledata.reads[1], l1=input_reads[0], l2=input_reads[1])
        srst2_log = self.run_quay_image(_TOOL, cmd=srst2_cmd_full)
        print(srst2_log)
        sampledata.reference_nfasta = mlst_db_abs
        sampledata.mlst_results = "{}__mlst__{}_{}__results.txt".format(sampledata.prefix, genus, species)
    # Orthologs-based phylogenetic tree construction
    def run_orthomcl(self):
        pass
    def handle(self, sdarr: SampleDataArray):
        functions = (self.run_fastqc, self.run_trimmomatic, self.run_cutadapt, self.run_spades, self.run_prokka,
                     self.run_srst2)
        for idx, func in enumerate(functions):
            logging.info("Starting the pipeline step {} of {}".format(idx, len(functions)))
            _ = Utils.single_core_queue(func, sdarr.lines)


class Utils:
    @staticmethod
    def single_core_queue(func, queue: list):
        return [func(i) for i in queue]
    @staticmethod
    def multi_core_queue(func, queue: list, processes: int = int(subprocess.getoutput("nproc").strip())):
        import multiprocessing
        pool = multiprocessing.Pool(processes=processes)
        output = pool.map(func, queue)
        pool.close()
        pool.join()
        return output


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    validator = ArgValidator()
    sampleDataArray = SampleDataArray.parse(validator.sampledata_file)
    handler = Handler(validator.output_dir)
    logging.INFO("The pipeline processing has been started")
    handler.handle(sampleDataArray)
    logging.INFO("The pipeline processing has been completed")
