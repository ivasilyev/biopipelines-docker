#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOT to be launched with Docker
# Get command:
# curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py

import os
import re
import logging
import subprocess
import multiprocessing
from time import sleep


class ArgValidator:
    def __init__(self):
        import argparse
        _STEPS = list(range(1, len(Handler.TOOLS) + 1))
        _STAGES = "<{}>".format("|".join([str(i) for i in _STEPS]))
        _STAGES_DESCRIPTION = "\n".join(["{}. {};".format(idx + 1, i) for idx, i in enumerate(Handler.TOOLS)])
        parser = argparse.ArgumentParser(description="""
Run prokaryotic genome analysis pipeline for group of files with given taxa information""".strip(),
                                         epilog="""
Stages: https://github.com/boulygina/bioinformatics-pipelines/blob/master/Prokaryotes_analysis/prokaryotes.md
Description:
{}
            """.format(_STAGES_DESCRIPTION).strip())
        parser.add_argument(
            "-i", '--input', metavar='<input.sampledata>', required=True,
            help="""
A tab-delimited table file containing information for each sample (strain) per row. 
Columns: 
1. Short sample name or strain, e.g. 'sample_01' or 'eco_O157:H7'; 
2. Paired end reads divided with semicolon, e.g. '/reads/eco_01_R1.fastq.gz;/reads/eco_01_R2.fastq.gz'; 
3. Taxon information divided with spaces, e.g. 'Escherichia coli O157:H7'""".strip())
        parser.add_argument("-s", "--start", help="Stage to start the pipeline, inclusive", type=int, default=_STEPS[0],
                            metavar=_STAGES, choices=_STEPS)
        parser.add_argument('-f', '--finish', help='Stage to finish the pipeline, inclusive', type=int,
                            default=_STEPS[-1], metavar=_STAGES, choices=_STEPS)
        parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Output directory', required=True)
        self._namespace = parser.parse_args()
        self.sampledata_file = self._namespace.input
        self.threads = multiprocessing.cpu_count()
        self.stages_to_do = []
        self.output_dir = self._namespace.output_dir
        self.log_dir = os.path.join(self.output_dir, "log", Utils.get_time())

    def validate(self):
        start_point = self._namespace.start
        finish_point = self._namespace.finish
        if start_point > finish_point:
            Utils.log_and_raise("The start stage number ({}) must be less than the end stage number ({})".format(
                start_point, finish_point))
        # Make 'stages_to_do' zero-based
        self.stages_to_do = range(start_point - 1, finish_point)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        else:
            logging.warning("The output directory exists: '{}'".format(self.output_dir))


class SampleDataArray:
    lines = []

    def validate(self):
        self.lines = [i for i in self.lines if i.exists]

    @staticmethod
    def parse(file):
        arr = SampleDataArray()
        with open(file, mode="r", encoding="utf-8") as f:
            arr.lines = [SampleDataLine.parse(j) for j in [i.strip() for i in f] if len(j) > 0]
            arr.validate()
            f.close()
        return arr


class SampleDataLine:
    exists = False
    prefix, genome, plasmid, genome_genbank, genome_pfasta, reference_nfasta, mlst_results = ("",) * 7

    def __init__(self, sample_name: str, sample_reads: list, taxa: list):
        # e.g "ecoli_sample", ["reads.1.fq", "reads.2.fq"], ["Escherichia", "coli", "O157:H7"]
        self.name = sample_name.strip()
        self.reads = [i.strip() for i in sample_reads]
        if not all([os.path.isfile(i) for i in self.reads]):
            logging.warning("File(s) not found: '{}'".format(self.reads))
        else:
            self.exists = True
            self.extension = self.get_extension(self.reads[0].strip())
            self.taxa = self._parse_taxa(taxa)
            self._set_prefix()

    @staticmethod
    def parse(line: str):
        name, reads, taxa = [j for j in [i.strip() for i in line.split("\t")] if len(j) > 0]
        reads = SampleDataLine._parse_reads(reads.split(";"))
        if len(reads) > 2:
            logging.warning("More than 2 paired end read files were given, only 2 first will be used: '{}'".format(
                reads))
            reads = reads[:2]
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
        return sorted([j for j in [i.strip() for i in reads] if len(j) > 0])

    @staticmethod
    def _parse_taxa(taxa: list):
        out = {i: "" for i in ("genus", "species", "strain")}
        taxa = [j for j in [i.strip() for i in taxa] if len(j) > 0]
        if len(taxa) == 0:
            Utils.log_and_raise("Empty taxon data!")
        out["genus"] = taxa[0].capitalize()
        if len(taxa) > 1:
            if not any(i.isdigit() for i in taxa[1]):
                sp = taxa[1].replace(".", "").lower()
                if sp != "sp":
                    out["species"] = sp
            else:
                out["strain"] = taxa[1]
        if len(taxa) > 2:
            out["strain"] = taxa[2]
        return out

    def _set_prefix(self):
        _MAX_PREFIX_LENGTH = 4
        if len(self.taxa["species"]) >= _MAX_PREFIX_LENGTH - 1:
            self.prefix = self.taxa["genus"][0].upper() + self.taxa["species"][:_MAX_PREFIX_LENGTH - 1].lower()
        else:
            self.prefix = self.taxa["genus"][:_MAX_PREFIX_LENGTH - len(
                self.taxa["species"])].capitalize() + self.taxa["species"].lower()

    def set_reads(self, reads: list):
        reads = self._parse_reads(reads)
        self.reads = reads


class Handler:
    TOOLS = ("fastqc", "trimmomatic", "cutadapt", "spades", "prokka", "bowtie2", "samtools", "vcftools", "snpeff",
             "srst2", "extract_pfasta_from_gbk", "orthomcl")

    def __init__(self, output_dir: str):
        self.output_dir_root = os.path.normpath(output_dir)
        if os.path.exists(self.output_dir_root):
            logging.warning("The path exists: '{}'".format(self.output_dir_root))
        os.makedirs(self.output_dir_root, exist_ok=True)
        # Output paths for each step
        self.output_dirs = {i: os.path.normpath(os.path.join(self.output_dir_root, "{}_{}".format(idx + 1, i))) for
                            idx, i in enumerate(self.TOOLS)}
        self.output_dirs = {i: os.path.normpath(os.path.join(
            self.output_dir_root, "_".join([str(idx + 1).zfill(len(str(len(self.TOOLS)))), i]))) for
            idx, i in enumerate(self.TOOLS)}
        # TODO Implenent these lines instead of method constants
        self._current_tool_name = ""
        self._current_tool_attempt = 0
        self._current_tool_max_attempts = 5

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

    def run_quay_image(self, img_name, img_tag: str = None, repo_name: str = "biocontainers", cmd: str = "echo",
                       bad_phrases: list = (), attempts: int = 5):
        import json
        # Get API response
        if not img_tag:
            attempt = 0
            url = "https://quay.io/api/v1/repository/{}/{}".format(repo_name, img_name)
            while attempt <= attempts:
                attempt += 1
                try:
                    api_response = json.loads(self.get_page(url))
                    img_tag = sorted(set(api_response.get("tags")))[-1]
                    break
                except json.decoder.JSONDecodeError:
                    logging.warning("Cannot get API response from the URL '{}' for attempt {} of {}".format(
                        url, attempt, attempts))
            if attempt > attempts:
                logging.warning("Exceeded attempts number to get API response from the URL '{}'".format(url))
        # Pull & run image
        img_name_full = "quay.io/{}/{}:{}".format(repo_name, img_name, img_tag)
        return Utils.run_image(img_name=img_name_full, container_cmd=cmd, bad_phrases=bad_phrases, attempts=attempts)

    @staticmethod
    def clean_path(path):
        if os.path.exists(path):
            logging.warning("The path exists and will be replaced with the new data: '{}'".format(path))
            if path != "/":  # I know this is useless
                logging.debug(subprocess.getoutput("echo Remove '{o}'; rm -rf {o}".format(o=os.path.normpath(path))))
        os.makedirs(path, exist_ok=True)  # path could be created under root and could not be deleted easily

    @staticmethod
    def process_reads(sampledata: SampleDataLine, out_dir: str, suffix: str):
        out_dir, suffix = [i.strip() for i in (out_dir, suffix)]
        processed_reads = [
            os.path.join(out_dir, "{}_{}.{}{}".format(sampledata.name, suffix, idx + 1, sampledata.extension)) for
            idx, i in enumerate(sampledata.reads)]
        return processed_reads

    # Pipeline steps
    def run_fastqc(self, sampledata: SampleDataLine, skip: bool = False):
        _TOOL = "fastqc"
        # One per read sample
        for idx, reads_file in enumerate(sampledata.reads):
            stage_dir = os.path.join(self.output_dirs[_TOOL], sampledata.name,
                                     "{}_{}".format(sampledata.name, idx + 1))
            cmd = """
            bash -c \
                'cd {o};
                 {T} -t {c} {r} -o {o};
                 chmod -R 777 {o}'
            """.format(T=_TOOL, c=validator.threads, r=reads_file, o=stage_dir)
            if not skip:
                self.clean_path(stage_dir)
                log = self.run_quay_image(_TOOL, cmd=cmd)
                Utils.append_log(log, _TOOL, sampledata.name)
            else:
                logging.info("Skip.")
        # Reads are unchanged, so there is nothing to return

    def run_trimmomatic(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "trimmomatic"
        stage_dir = os.path.join(self.output_dirs[_TOOL], sampledata.name)
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        untrimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix="{}_untrinned".format(_TOOL))
        cmd = """
        bash -c \
            'cd {o};
             {T} PE -threads {t} -phred33 {r1} {r2} {t1} {u1} {t2} {u2} \
                ILLUMINACLIP:/data2/bio/ecoli_komfi/adapters.fasta:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
             chmod -R 777 {o}'
        """.format(T=_TOOL, t=validator.threads, o=stage_dir, r1=sampledata.reads[0], r2=sampledata.reads[1],
                   t1=trimmed_reads[0], t2=trimmed_reads[1], u1=untrimmed_reads[0], u2=untrimmed_reads[1])
        if not skip:
            self.clean_path(stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip.")
        sampledata.set_reads(trimmed_reads)

    def run_cutadapt(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "cutadapt"
        _ADAPTER = "AGATCGGAAGAG"
        stage_dir = os.path.join(self.output_dirs["cutadapt"], sampledata.name)
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        cmd = """
        bash -c \
            'cd o;
             {T} -a {a} -A {a} -m 50 -o {t1} -p {t2} {r1} {r2};
             chmod -R 777 {o}'
        """.format(T=_TOOL, a=_ADAPTER, r1=sampledata.reads[0], r2=sampledata.reads[1], t1=trimmed_reads[0],
                   t2=trimmed_reads[1], o=stage_dir)
        if not skip:
            self.clean_path(stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip.")
        sampledata.set_reads(trimmed_reads)

    def run_spades(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "spades"
        """
        # Sample launch:
        IMG=quay.io/biocontainers/spades:3.13.1--0 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG bash -c \
            'TOOL=$(find /usr/local/share/ -name spades.py | grep spades | head -n 1) && $TOOL -v'
        """
        stage_dir = os.path.join(self.output_dirs[_TOOL], sampledata.name)
        assemblies = {"genome": "", "plasmid": ""}
        for assembly_type in assemblies:
            assembly_dir = os.path.join(stage_dir, assembly_type)
            cmd_append = ""
            if assembly_type == "plasmid":
                cmd_append = " --plasmid"
            cmd = """
            bash -c \
                'cd {o};
                 TOOL=$(find /usr/local/share/ -name {t}.py | grep {t} | head -n 1) && \
                 $TOOL --careful -o {o} -1 {r1} -2 {r2}{a};
                 chmod -R 777 {o}'
            """.format(o=assembly_dir, t=_TOOL, r1=sampledata.reads[0], r2=sampledata.reads[1], a=cmd_append)
            if not skip:
                self.clean_path(assembly_dir)
                log = self.run_quay_image(_TOOL, cmd=cmd)
                Utils.append_log(log, _TOOL, sampledata.name)
            else:
                logging.info("Skip.")
            assemblies[assembly_type] = os.path.join(assembly_dir, "contigs.fasta")
        sampledata.genome = assemblies["genome"]
        sampledata.plasmid = assemblies["plasmid"]

    def run_prokka(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "prokka"
        """
        # Sample launch:
        IMG=quay.io/biocontainers/prokka:1.12--pl526_0 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG prokka
        """
        stage_dir = os.path.join(self.output_dirs[_TOOL], sampledata.name)
        taxa_append = ""
        for taxon_name in ("genus", "species", "strain"):
            taxon_value = sampledata.taxa.get(taxon_name)
            if len(taxon_value) > 0:
                taxa_append = "{} --{} {}".format(taxa_append, taxon_name, taxon_value)
        cmd = """
        bash -c \
            'cd {o}
             {T} --compliant --centre UoN --cpu {c} --outdir {o} --force --prefix {n} --locustag {n} {a} {g}
             chmod -R 777 {o}'
        """.format(T=_TOOL, g=sampledata.genome, a=taxa_append, n=sampledata.name, o=stage_dir,
                   c=validator.threads)
        if not skip:
            self.clean_path(stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip.")
        sampledata.genome_genbank = os.path.join(stage_dir, "{}.gbk".format(sampledata.name))
        sampledata.genome_pfasta = os.path.join(stage_dir, "{}.faa".format(sampledata.name))

    # SNP calling
    def run_bowtie2(self, sampledata: SampleDataLine, skip: bool = False):
        pass

    def run_samtools(self, sampledata: SampleDataLine, skip: bool = False):
        pass

    def run_vcftools(self, sampledata: SampleDataLine, skip: bool = False):
        pass

    # SNP annotation
    def run_snpeff(self, sampledata: SampleDataLine, skip: bool = False):
        pass

    # MLST typing
    def run_srst2(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample, full cleaning is NOT required
        _TOOL = "srst2"
        _GETMLST_ATTEMPTS = 5
        _SRST2_ATTEMPTS = 5
        """
        # Sample launch:
        IMG=quay.io/biocontainers/srst2:0.2.0--py27_2 && \ 
        docker pull $IMG && \
        docker run --rm --net=host -it $IMG getmlst.py -h
        
        docker run --rm --net=host -it $IMG srst2 -h
        """
        tool_dir = self.output_dirs[_TOOL]
        stage_dir = os.path.join(tool_dir, sampledata.name)
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

        srst2 --output test --input_pe *.fastq.gz --mlst_db Klebsiella_pneumoniae.fasta \
            --mlst_definitions kpneumoniae.txt --mlst_delimiter '_'
        
        # Sample srst2 output:
        
        Klebsiella_MLST__mlst__Klebsiella_pneumoniae__results.txt
        """
        # Default cmd
        srst2_cmd = """
        srst2 --output test --input_pe *.fastq.gz --mlst_db {} --mlst_definitions {} --mlst_delimiter '_'
        """.format(mlst_db_local, mlst_definitions_local)
        sampledata.reference_nfasta = mlst_db_abs
        sampledata.mlst_results = "{}__mlst__{}_{}__results.txt".format(sampledata.name, genus, species)
        if skip:
            logging.info("Skip.")
            return
        if not all(os.path.isfile(i) for i in (mlst_db_abs, mlst_definitions_abs)):
            os.makedirs(tool_dir, exist_ok=True)
            getmlst_attempt = 0
            while getmlst_attempt < _GETMLST_ATTEMPTS:
                getmlst_attempt += 1
                getmlst_cmd = """
                bash -c \
                    'cd {o};
                     getmlst.py --species "{g} {s}";
                     chmod -R 777 {o}'
                """.format(g=genus, s=species, o=tool_dir)
                getmlst_log = self.run_quay_image("srst2", cmd=getmlst_cmd)
                srst2_cmd = getmlst_log.split("Suggested srst2 command for use with this MLST database:")[-1].strip()
                Utils.append_log(getmlst_log, "getmlst", sampledata.name)
                if not srst2_cmd.startswith("srst2"):
                    logging.warning("`getmlst.py` did not finish correctly for attempt {} of {}".format(
                        getmlst_attempt, _GETMLST_ATTEMPTS))
                else:
                    break
            if getmlst_attempt == _GETMLST_ATTEMPTS:
                logging.warning("Exceeded attempts number for `getmlst.py` to finish processing correctly")
        # The input read files must be named by a strict pattern:
        # https://github.com/katholt/srst2#input-read-formats-and-options
        input_reads = ["{}_{}.fastq.gz".format(sampledata.name, idx + 1) for idx, i in enumerate(sampledata.reads)]
        srst2_cmd = srst2_cmd.replace(
            "*.fastq.gz", " ".join(input_reads)).replace(
            "--output test", "--output {}_MLST".format(sampledata.name)).replace(
            mlst_db_local, mlst_db_abs).replace(
            mlst_definitions_local, mlst_definitions_abs)
        srst2_cmd_full = """
        bash -c \
            'cd {o};
             ln -s {r1} {l1};
             ln -s {r2} {l2};
             {c} --log --threads {t};
             chmod -R 777 {o}'
        """.format(c=srst2_cmd.strip(), o=stage_dir, t=validator.threads,
                   r1=sampledata.reads[0], r2=sampledata.reads[1], l1=input_reads[0], l2=input_reads[1])
        self.clean_path(stage_dir)
        srst2_log = self.run_quay_image(_TOOL, cmd=srst2_cmd_full, attempts=_SRST2_ATTEMPTS,
                                        bad_phrases=["Encountered internal Bowtie 2 exception",
                                                     "[main_samview] truncated file."])
        Utils.append_log(srst2_log, _TOOL, sampledata.name)

    def extract_pfasta_from_gbk(self, sampledata: SampleDataLine, skip: bool = False):
        _TOOL = "extract_pfasta_from_gbk"
        tool_dir = self.output_dirs[_TOOL]
        os.makedirs(tool_dir, exist_ok=True)
        sampledata.genome_pfasta = os.path.join(tool_dir, "{}.genome.protein.fasta".format(sampledata.name))
        if not skip:
            log = Utils.run_image("ivasilyev/orthomcl-mysql:latest", container_cmd="""
            python3 /opt/my_tools/{}.py -i {} -s {} -o {}
            """.format(_TOOL, sampledata.genome_genbank, sampledata.name, sampledata.genome_pfasta))
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip.")

    # Orthologs-based phylogenetic tree construction
    def run_orthomcl(self, sampledata_array: SampleDataArray, skip: bool = False):
        # One per all samples
        _TOOL = "orthomcl"
        tool_dir = self.output_dirs[_TOOL]
        if skip:
            logging.info("Skip.")
            return
        self.clean_path(tool_dir)
        orthomcl_sampledata = os.path.join(tool_dir, "orthomcl_abbreviations.sampledata")
        with open(orthomcl_sampledata, mode="w", encoding="utf-8") as f:
            for sampledata in sampledata_array.lines:
                f.write("{}\t{}\n".format(sampledata.genome_pfasta, sampledata.name))
            f.close()
        log = Utils.run_image(img_name="ivasilyev/orthomcl-mysql:latest",
                              container_cmd="""
                              bash -c \
                                'service mysql restart;
                                 python3 /opt/my_tools/pipeline_handler.py -i {s} -o {o};
                                 chmod -R 777 {o}'
                              """.format(s=orthomcl_sampledata, o=tool_dir))
        Utils.append_log(log, _TOOL, "all")

    def handle(self, sampledata_array: SampleDataArray):
        _SAMPLE_METHODS = (self.run_fastqc, self.run_trimmomatic, self.run_cutadapt, self.run_spades, self.run_prokka,
                           self.run_bowtie2, self.run_samtools, self.run_vcftools, self.run_snpeff, self.run_srst2,
                           self.extract_pfasta_from_gbk)
        _GROUP_METHODS = (self.run_orthomcl,)
        for idx, func in enumerate(_SAMPLE_METHODS + _GROUP_METHODS):
            try:
                logging.info("Starting the pipeline step {} of {} ({} in total)".format(
                    idx + 1, len(_SAMPLE_METHODS + _GROUP_METHODS), len(validator.stages_to_do)))
                # Per-sample processing
                if func in _SAMPLE_METHODS:
                    queue = [(func, i, idx not in validator.stages_to_do) for i in sampledata_array.lines]
                    _ = Utils.single_core_queue(Utils.wrap_func, queue)
                # Per-group functions
                else:
                    func(sampledata_array, skip=idx not in validator.stages_to_do)
            except PermissionError:
                logging.critical("Cannot process the step {}, please run the command 'sudo chmod -R 777 {}'".format(
                    idx, self.output_dir_root))


class Utils:
    @staticmethod
    def get_time():
        from datetime import datetime
        now = datetime.now()
        output_list = []
        for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
            time_unit = str(time_unit)
            if len(time_unit) < 2:
                time_unit = '0' + time_unit
            output_list.append(time_unit)
        return '-'.join(output_list)

    @staticmethod
    def remove_empty_values(input_list):
        return [j for j in [i.strip() for i in input_list] if len(j) > 0]

    @staticmethod
    def log_and_raise(msg):
        logging.critical(msg)
        raise ValueError(msg)

    @staticmethod
    def append_log(msg: str, tool_name: str, sample_name: str):
        file = os.path.join(validator.log_dir, "{}_{}.log".format(tool_name, sample_name))
        with open(file, mode="a", encoding="utf-8") as f:
            f.write(msg + "\n")
            f.close()

    @staticmethod
    def is_log_valid(log: str, bad_phrases: list):
        log_lines = Utils.remove_empty_values(log.replace("\r", "\n").split("\n"))
        bad_phrases = Utils.remove_empty_values(bad_phrases)
        if len(bad_phrases) == 0:
            return True
        return all([i not in j for i in bad_phrases for j in log_lines])

    @staticmethod
    def run_until_valid_output(cmd: str, bad_phrases: list, attempts: int = 5, ping_required: bool = False):
        attempt = 0
        log = ""
        while attempt < attempts:
            attempt += 1
            log_cmd = re.sub("[\r\n ]+", " ", cmd.strip())
            logging.debug("Executing the command: `{}`".format(log_cmd))
            log = subprocess.getoutput(cmd)
            if not Utils.is_log_valid(log, bad_phrases):
                logging.warning("An error phrase was found in log for attempt {} of {}.".format(attempt, attempts))
                sleep(5)
                if ping_required:
                    # Ping Google just to keep the node DNS working
                    _ = subprocess.getoutput("ping -c 10 google.com")
            else:
                break
            if attempt == attempts:
                logging.warning("Exceeded attempts number to get execution output without failure messages. "
                                "The command seems to be not finished correctly: `{}`".format(log_cmd))
        return log

    @staticmethod
    def run_image(img_name: str, container_cmd: str, bad_phrases: list = (), attempts: int = 5):
        _DOCKER_RUN_CMD = "docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it"
        _COMMON_PHRASES = ["Error response from daemon", ]
        logging.info("Using image: '{}'".format(img_name))
        bad_phrases = list(bad_phrases) + _COMMON_PHRASES
        docker_cmd = "docker pull {a} && {b} {a}".format(a=img_name, b=_DOCKER_RUN_CMD)
        # cmd may contain curly braces, so str.format() is not usable
        out_cmd = docker_cmd.strip() + " " + container_cmd.strip()
        return Utils.run_until_valid_output(cmd=out_cmd, bad_phrases=bad_phrases, attempts=attempts, ping_required=True)

    @staticmethod
    def wrap_func(args: list):
        return args[0](*args[1:])

    @staticmethod
    def single_core_queue(func, queue: list):
        return [func(i) for i in queue]

    @staticmethod
    def multi_core_queue(func, queue: list, processes: int = None):
        if not processes:
            processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=processes)
        output = pool.map(func, queue)
        pool.close()
        pool.join()
        return output


if __name__ == '__main__':
    validator = ArgValidator()
    mainLogFile = os.path.join(validator.log_dir, "main.log")
    os.makedirs(validator.log_dir, exist_ok=True)
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s",
                        handlers=[logging.FileHandler(mainLogFile), logging.StreamHandler()])
    validator.validate()
    sampleDataArray = SampleDataArray.parse(validator.sampledata_file)
    handler = Handler(validator.output_dir)
    logging.info("The pipeline processing has been started")
    handler.handle(sampleDataArray)
    logging.info("The pipeline processing has been completed")
