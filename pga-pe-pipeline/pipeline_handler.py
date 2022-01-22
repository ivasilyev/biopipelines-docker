#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOT to be launched with Docker
# Get command:
# curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py

import os
import re
import json
import logging
import subprocess
import multiprocessing
from shutil import copy2
from datetime import datetime
from itertools import product
from collections.abc import Mapping
from time import sleep, perf_counter


BOWTIE2_HG_IDX_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz"
PULL_RETRIES = 5


class ArgValidator:
    def __init__(self):
        import argparse
        _handler = Handler()
        _handler_methods = [i.__name__ for i in (_handler.sample_methods + _handler.group_methods)]
        _STEPS = list(range(1, len(_handler_methods) + 1))
        _STAGES = "<{}>".format("|".join([str(i) for i in _STEPS]))
        _STAGES_DESCRIPTION = "\n".join(["{}. {};".format(idx + 1, i) for idx, i in enumerate(_handler_methods)])

        parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                         description="""
Run prokaryotic genome analysis pipeline for group of raw read files with given taxa information""".strip(),
                                         epilog="""
Stages: https://github.com/boulygina/bioinformatics-pipelines/blob/master/Prokaryotes_analysis/prokaryotes.md

Description:
{}
            """.format(_STAGES_DESCRIPTION).strip())
        parser.add_argument(
            "-i", "--input", metavar="<input.sampledata>", required=True,
            help="""
A tab-delimited table file containing information for each sample (strain) per row. 
Columns: 
1. Short sample name or strain, e.g. 'sample_01' or 'eco_O157:H7'; 
2. Paired end reads divided with semicolon, e.g. '/reads/eco_01_R1.fastq.gz;/reads/eco_01_R2.fastq.gz'; 
3. Taxon information divided with spaces, e.g. 'Escherichia coli O157:H7'""".strip())
        parser.add_argument("-s", "--start", type=int, default=_STEPS[0], metavar=_STAGES,
                            choices=_STEPS, help="(Optional) Stage to start the pipeline, inclusive")
        parser.add_argument("-f", "--finish", type=int, default=_STEPS[-1], metavar=_STAGES,
                            choices=_STEPS, help="(Optional) Stage to finish the pipeline, inclusive")
        parser.add_argument("--hg_dir", metavar="<dir>", default="",
                            help="(Optional) Directory containing human genome bowtie2 index ('*.bt2') required for human genome decontamination")
        parser.add_argument("--blast_references", metavar="<int>", type=int, default=100,
                            help="(Optional) Number of BLAST references to fetch")
        parser.add_argument("--card_json", metavar="<card.json>", default="",
                            help="(Optional) CARD reference JSON required by the RGI")
        parser.add_argument("--refdata", metavar="<file>", default=(), nargs="+",
                            help="(Optional) Path(s) to RefData JSONs made by the 'cook_the_reference.py' required for alignment and coverage extrction")
        parser.add_argument("-o", "--output_dir", metavar="<dir>", required=True,
                            help="Output directory")

        self._namespace = parser.parse_args()
        self.sampledata_file = os.path.realpath(self._namespace.input)

        self.threads = multiprocessing.cpu_count()
        self.stages_to_do = []

        self.hg_index_dir = self._namespace.hg_dir
        if len(self.hg_index_dir) > 0:
            self.hg_index_dir = os.path.realpath(self.hg_index_dir)

        self.blast_reference_number = self._namespace.blast_references
        self.card_json = self._namespace.card_json

        self.refdata_files = Utils.remove_empty_values(self._namespace.refdata)
        if len(self.refdata_files) > 0:
            self.refdata_files = tuple([os.path.realpath(i) for i in self.refdata_files])

        self.output_dir = os.path.realpath(self._namespace.output_dir)
        self.log_dir = os.path.join(self.output_dir, "logs", Utils.get_time())

    def validate(self):
        start_point = self._namespace.start
        finish_point = self._namespace.finish
        if start_point > finish_point:
            Utils.log_and_raise("The start stage number ({}) must be less than the end stage number ({})".format(
                start_point, finish_point))
        # Make 'stages_to_do' zero-based
        self.stages_to_do = range(start_point - 1, finish_point)


class SampleDataLine:
    def __init__(self, sample_name: str, sample_reads: list, taxa):
        self.state = dict()
        self.prefix = ""
        self.chromosome_assembly = ""
        self.plasmid_assembly = ""
        self.genome_assembly = ""
        self.genbank = ""
        self.genomes = dict()
        self.faa = ""
        self.chromosome_annotation = ""
        self.blast_result_json = ""
        self.blast_result_table = ""

        self.reference_fna = ""
        self.srst2_result_table = ""

        # E.g "ecoli_sample", ["reads.1.fq", "reads.2.fq"], "Escherichia coli O157:H7"]
        self.name = sample_name.strip()

        self.reads = []
        self.raw_reads = []
        self.set_reads(sample_reads)
        self.raw_reads = list(self.reads)

        self.is_valid = False
        self._validate()
        self.extension = Utils.get_reads_file_extension(self.raw_reads[0])
        self.taxa_genus, self.taxa_species, self.taxa_strain = ["", ] * 3
        self.closest_reference_genbank = ""
        self._parse_taxa(taxa)

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def set_reads(self, reads: list):
        self.reads = sorted(Utils.remove_empty_values(reads))
        if len(self.reads) > 2:
            logging.warning("Too much of reads are given (the maximum is 2)!")
            self.reads = self.reads[:2]

    def _validate(self):
        c = 0
        # Validate reads
        if len(self.reads) < 2:  # PE are the only library strategy supported
            logging.warning("Not enough of reads are given (the minimum is 2)!")
            c += 1
        if not all([Utils.is_file_valid(i) for i in self.reads]):
            logging.warning("Some read files are missing!")
            c += 1
        self.is_valid = c == 0

    @staticmethod
    def parse(d: dict):
        """
        :param d: dict
        {"key_1": "sample_name", "key_2": "raw_reads", "key_1": "taxa", ...}
        """
        _KEYS = ("sample_name", "sample_reads", "taxa")
        if all(i in d.keys() for i in _KEYS):
            return SampleDataLine(d["sample_name"], d["sample_reads"], d["taxa"])
        keys = list(d.keys())
        return SampleDataLine(*[d[i] for i in keys[:3]])

    @property
    def reads_string(self):
        return " ".join(self.reads)

    @property
    def taxa_string(self):
        return " ".join([self.taxa_genus, self.taxa_species, self.taxa_strain])

    @property
    def is_taxa_valid(self):
        return len(self.taxa_string) > 0

    def set_taxa(self, genus: str = "", species: str = "", strain: str = ""):
        self.taxa_genus, self.taxa_species, self.taxa_strain = genus, species, strain
        if not self.is_taxa_valid:
            logging.warning("No taxonomy data specified, some steps will be passed for the sample '{}'".format(self.name))
        else:
            self._set_prefix()

    def _parse_taxa(self, taxa):
        self.set_taxa(**Utils.parse_taxa(taxa))

    def _set_prefix(self):
        _MAX_PREFIX_LENGTH = 4
        if not self.is_taxa_valid:
            return
        if len(self.taxa_species) >= _MAX_PREFIX_LENGTH - 1:
            self.prefix = "{}{}".format(self.taxa_genus[0].upper(),
                                        self.taxa_species[:_MAX_PREFIX_LENGTH - 1].lower())
        else:
            self.prefix = "{}{}".format(self.taxa_genus[:_MAX_PREFIX_LENGTH - len(self.taxa_species)].capitalize(),
                                        self.taxa_species.lower())


class SampleDataArray:
    def __init__(self):
        self.lines = dict()
        self.srst2_merged_table = ""
        self.blast_merged_table = ""
        self.roary_edited_newick = ""

    @property
    def blast_result_tables(self):
        return Utils.remove_empty_values([i.blast_result_table for i in self.lines.values()])

    @property
    def srst2_result_tables(self):
        return Utils.remove_empty_values([i.srst2_result_table for i in self.lines.values()])

    def validate(self):
        d = dict()
        for key in self.lines:
            line = self.lines[key]
            if line.is_valid:
                d[key] = line
            else:
                logging.warning("Invalid sample data for the sample: '{}'".format(line.name))
        if len(d.keys()) == 0:
            Utils.log_and_raise("No valid sample data lines, exit")
        self.lines = d

    def __len__(self):
        return len(self.lines)

    @staticmethod
    def parse(d: dict):
        """
        :param d: dict
        {sample_name: , sample_reads: list, taxa)}
        :return:
        """

        arr = SampleDataArray()
        arr.lines = {k: SampleDataLine.parse(v) for k, v in d.items()}
        arr.validate()
        return arr

    @staticmethod
    def load(file: str):
        """
        :param file: str

        An example sample data JSON:
        { "ecoli_sample":
            { "name": "ecoli_sample",
              "raw_reads": [ "reads.1.fq", "reads.2.fq" ],
              "taxa": "Escherichia coli O157:H7" }, }

        or

        { "ecoli_sample":
            { "name": "ecoli_sample",
              "raw_reads": [ "reads.1.fq", "reads.2.fq" ],
              "taxa": { "genus": "Escherichia",
                        "species": "coli",
                        "strain": "O157:H7", }, }, }

        :return: dict
        """
        return SampleDataArray.parse(Utils.load_dict(file))

    def export(self):
        return {k: v.reads for k, v in self.lines.items()}

    def to_2d_array(self):
        """
        :return: 2d array
        [
            [ "ecoli_sample", "reads.1.fq", "reads.2.fq" ],
            ...
        ]
        """
        return [[k] + v.reads for k, v in self.lines.items()]

    def dump(self, file: str):
        Utils.dump_dict(self.export(), file)


class Handler:
    def __init__(self, output_dir: str = ""):
        self.sample_methods = [
            self.run_fastqc,
            self.run_trimmomatic,
            self.run_cutadapt,
            self.remove_hg,
            self.run_spades,
            self.run_plasmid_merger,
            self.run_blast,
            self.run_quast,
            self.run_prokka,
            self.run_mgefinder,
            self.run_rgi,
            self.run_srst2,
        ]
        self.group_methods = [
            self.merge_srst2_results,
            self.merge_blast_results,
            self.run_roary,
            self.run_nbee_with_annotation,
        ]
        #
        self.valid = False
        self.output_dir_root = output_dir.strip()
        self.output_dirs = dict()
        #
        self._reference_dir = os.path.join(self.output_dir_root, "references")
        self.human_genome_reference_dir = os.path.join(self._reference_dir, "human_genome")
        self.blast_reference_dir = os.path.join(self._reference_dir, "blast")
        self.srst2_reference_dir = os.path.join(self._reference_dir, "srst2")
        self.mgefinder_reference_dir = os.path.join(self._reference_dir, "mgefinder")
        #
        self.card_reference_dir = os.path.join(self._reference_dir, "card")
        self.card_reference_json = ""
        #
        self.roary_reference_dir = os.path.join(self._reference_dir, "roary")
        #
        self._state = dict()
        #
        self.prepare_environment()

    def prepare_environment(self):
        if len(self.output_dir_root) == 0:
            return
        self.output_dir_root = os.path.normpath(self.output_dir_root)
        if os.path.exists(self.output_dir_root):
            logging.warning("The path exists: '{}'".format(self.output_dir_root))
        os.makedirs(self.output_dir_root, exist_ok=True)
        # Output paths for each step
        methods = self.sample_methods + self.group_methods
        for idx, method in enumerate(methods):
            prefix = str(idx + 1).zfill(len(str(len(methods))))
            method_name = "_".join([i for i in method.__name__.split("_") if len(i) > 0 and i != "run"])
            dir_name = os.path.normpath(os.path.join(self.output_dir_root, "_".join([prefix, method_name])))
            self.output_dirs[method.__name__] = dir_name
        self.valid = True

    def get_latest_quay_tag(self, repo_name, img_name):
        url = "https://quay.io/api/v1/repository/{}/{}".format(repo_name, img_name)
        api_response = json.loads(Utils.get_page(url))
        tool_tags = api_response["tags"].values()
        # Sample datetime:
        # Wed, 11 Apr 2018 17:20:46 -0000
        # https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior
        _ = [i.update(
            {"datetime": datetime.strptime(i["last_modified"], "%a, %d %b %Y %H:%M:%S %z")})
             for i in tool_tags]
        img_tag = sorted(tool_tags, key=lambda x: x["datetime"], reverse=True)[0]["name"]
        return img_tag

    def run_quay_image(
            self, img_name, img_tag: str = None, repo_name: str = "biocontainers",
            cmd: str = "echo", sample_name: str = "all", bad_phrases: list = (), attempts: int = 5
    ):
        if not img_tag:
            # Get API response
            attempt = 0
            while attempt <= attempts:
                attempt += 1
                try:
                    img_tag = self.get_latest_quay_tag(repo_name, img_name)
                    break
                except json.decoder.JSONDecodeError:
                    logging.warning(
                        f"Cannot get API response for the image the URL '{img_name}' for attempt {attempt} of {attempts}"
                    )
            if attempt > attempts:
                logging.warning(
                    f"Exceeded attempts number to get API response for the image '{img_name}'"
                )
        # Pull & run image
        img_name_full = "quay.io/{}/{}:{}".format(repo_name, img_name, img_tag)

        # Update software
        self.update_state({
            sample_name: {
                "software": {
                    repo_name: {
                        img_name: img_name_full
                    }
                }
            }
        })
        return Utils.run_image(
            img_name=img_name_full, container_cmd=cmd, bad_phrases=bad_phrases, attempts=attempts
        )

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

    # Single processing pipeline steps

    def _run_fastqc(self, sampledata_name: str, reads_file, out_dir: str):
        # One per read file, two per sample
        _TOOL = "fastqc"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1 && \
        docker pull ${IMG} && \
        docker run --rm -v /data1/data1 --net=host -it ${IMG} bash
        """
        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd "{out_dir}";
            {_TOOL} \
                --extract \
                --outdir "{out_dir}" \
                --threads {argValidator.threads} \
                "{reads_file}";
            rm -f *.zip;
            chmod -R a+rw "{out_dir}";
        '
        """
        self.clean_path(out_dir)
        log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata_name)
        Utils.append_log(log, _TOOL, sampledata_name)

    @staticmethod
    def _parse_fastqc_result(directory: str):
        """
        :param directory: str
        :return: dict

        'summary.txt' contains 3 columns: the value column (PASS|WARN|FAIL), the key column.
        The 3rd column contains only the repeating FastQC input reads file name.

        The key column:

        Basic Statistics
        Per base sequence quality
        Per tile sequence quality
        Per sequence quality scores
        Per base sequence content
        Per sequence GC content
        Per base N content
        Sequence Length Distribution
        Sequence Duplication Levels
        Overrepresented sequences
        Adapter Content
        """
        summary_table = Utils.locate_file_by_tail(directory, "summary.txt")
        if len(summary_table) == 0:
            return dict()
        arr = Utils.load_2d_array(summary_table)
        return {i[1]: i[0] for i in arr if len(i) > 1}

    def run_fastqc(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        stage_name = Utils.get_caller_name()
        tool_dir = self.output_dirs[stage_name]
        stage_dir = os.path.join(tool_dir, sampledata.name)

        for index, reads_file in enumerate(sampledata.reads):
            sub_stage_dir = os.path.join(stage_dir, f"{sampledata.name}_{index + 1}")
            if not skip:
                self._run_fastqc(
                    sampledata_name=sampledata.name, reads_file=reads_file, out_dir=sub_stage_dir
                )
            fastqc_results = self._parse_fastqc_result(sub_stage_dir)
            if len(fastqc_results.keys()) == 0 and not skip:
                logging.warning("No FastQC results!")
            self.update_state({
                sampledata.name: {
                    stage_name: {
                        os.path.basename(reads_file): fastqc_results
                    }
                }
            })

    def run_trimmomatic(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "trimmomatic"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_dir = os.path.join(self.output_dirs[Utils.get_caller_name()], sampledata.name)
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix="{}_trimmed".format(_TOOL))
        untrimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix="{}_untrimmed".format(_TOOL))
        # There is no '-version' or '--help' CLI argument
        cmd = f"""
        bash -c '
            cd {stage_dir};
            {_TOOL} \
                PE \
                -threads {argValidator.threads} \
                -phred33 \
                {sampledata.reads_string} \
                {trimmed_reads[0]} {untrimmed_reads[0]} {trimmed_reads[1]} {untrimmed_reads[1]} \
                ILLUMINACLIP:/data2/bio/ecoli_komfi/adapters.fasta:2:30:10 \
                LEADING:3 \
                TRAILING:3 \
                SLIDINGWINDOW:4:15 \
                MINLEN:36;
             chmod -R a+rw {stage_dir};
        '
        """
        if not skip:
            self.clean_path(stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip {}".format(Utils.get_caller_name()))
        sampledata.set_reads(trimmed_reads)

    def run_cutadapt(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "cutadapt"
        _ADAPTER = "AGATCGGAAGAG"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/cutadapt:3.5--py37h73a75cf_0 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_dir = os.path.join(self.output_dirs[Utils.get_caller_name()], sampledata.name)
        trimmed_reads = self.process_reads(sampledata, out_dir=stage_dir, suffix=_TOOL)
        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd {stage_dir};
            {_TOOL} \
                --adapter {_ADAPTER} \
                -A {_ADAPTER} \
                --cores {argValidator.threads} \
                --minimum-length 50 \
                --output {trimmed_reads[0]} \
                --paired-output {trimmed_reads[1]} \
                {sampledata.reads_string};
            chmod -R a+rw {stage_dir};
        '
        """
        if not skip:
            self.clean_path(stage_dir)
            log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip {}".format(Utils.get_caller_name()))
        sampledata.set_reads(trimmed_reads)

    @staticmethod
    def _parse_bowtie2_index_mask(directory: str):
        indices = sorted(Utils.locate_file_by_tail(directory, ".bt2", True), key=len)
        if len(indices) > 0:
            return ".".join(indices[0].split(".")[:-2])
        return ""

    def _download_hg_reference(self, directory: str):
        self.clean_path(directory)
        cmd = f"""
        bash -c '
            cd {directory};
            curl -fsSL "{BOWTIE2_HG_IDX_URL}" | tar -xzf -;
        '
        """
        logging.info("Downloaded the reference human genome")
        return Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)

    @staticmethod
    def _parse_bowtie2_log(s: str):
        _KEYS = [
            " were paired", " aligned concordantly 0 times", " aligned concordantly exactly 1 time",
            " aligned concordantly >1 times",
            " pairs aligned concordantly 0 times", " aligned discordantly 1 time",
            " pairs aligned 0 times concordantly or discordantly", " mates make up the pairs",
            " aligned 0 times", " aligned exactly 1 time", " aligned >1 times",
            " overall alignment rate"
        ]
        lines = Utils.split_lines(s)
        return {i.strip(): re.sub(f"{i}.*$", "", j) for i, j in product(_KEYS, lines) if i in j}

    def remove_hg(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "bowtie2"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/bowtie2:2.4.4--py37h13ad519_0 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_name = Utils.get_caller_name()
        stage_dir = self.output_dirs[stage_name]
        state_file = os.path.join(stage_dir, "state.json")
        if Utils.is_file_valid(state_file):
            self.update_state({stage_name: Utils.load_dict(state_file)})
        if skip:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        if len(argValidator.hg_index_dir) > 0:
            self.human_genome_reference_dir = argValidator.hg_index_dir
        os.makedirs(self.human_genome_reference_dir, exist_ok=True)
        index_mask = self._parse_bowtie2_index_mask(self.human_genome_reference_dir)
        if len(index_mask) == 0:
            logging.warning("The directory does not contain valid bowtie2 (*.bt2) indexes: '{}'".format(self.human_genome_reference_dir))
            log = self._download_hg_reference(self.human_genome_reference_dir)
            Utils.append_log(log, _TOOL, sampledata.name)
            index_mask = self._parse_bowtie2_index_mask(self.human_genome_reference_dir)

        mapped_reads_dir = os.path.join(stage_dir, "mapped")
        mapped_reads_file = os.path.join(mapped_reads_dir, "{}.sam".format(sampledata.name))
        unmapped_reads_dir = os.path.join(stage_dir, "unmapped", sampledata.name)
        # bowtie2 may mess with double extensions, e.g. ".fq.1.gz" instead of ".1.fq.gz"
        unmapped_file_mask = os.path.join(unmapped_reads_dir, sampledata.name)
        os.makedirs(mapped_reads_dir, exist_ok=True)
        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd {stage_dir};
            {_TOOL} \
                --gbar 1 \
                --local \
                --mp 3 \
                --threads {argValidator.threads} \
                --un-conc {unmapped_file_mask} \
                -1 {sampledata.reads[0]} \
                -2 {sampledata.reads[1]} \
                -D 20 \
                -L 3 \
                -N 1 \
                -R 3 \
                -S {mapped_reads_file} \
                -x {index_mask};
            chmod -R a+rw {stage_dir}'
        """
        self.clean_path(unmapped_reads_dir)
        log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
        Utils.append_log(log, _TOOL, sampledata.name)
        unmapped_reads_files = [i for i in Utils.scan_whole_dir(unmapped_reads_dir) if unmapped_file_mask in i]
        if len(unmapped_reads_files) != 2:
            logging.warning(
                "The number of output unmapped reads is not equal to 2, please check the directory: '{}'".format(
                    unmapped_reads_dir))
        else:
            unmapped_reads_files_new = []
            # SPAdes needs clearly specified extension
            for unmapped_reads_file_old in unmapped_reads_files:
                unmapped_reads_file_new = "{}{}".format(unmapped_reads_file_old, sampledata.extension)
                os.rename(unmapped_reads_file_old, unmapped_reads_file_new)
                unmapped_reads_files_new.append(unmapped_reads_file_new)
            sampledata.set_reads(unmapped_reads_files_new)
        # Parse log
        state_dict = self._parse_bowtie2_log(log)
        Utils.dump_dict(state_dict, state_file)
        self.update_state({
            sampledata.name: {
                stage_name: state_dict
            }
        })

    def run_spades(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "spades"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/spades:3.13.1--0 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        
        # Tool executable lookup:
        export TOOL="$(find /usr/local/ -name "spades.py" -type f 2>/dev/null | grep 'spades.py$' | head -n 1)"
        echo "${TOOL}"
        """
        stage_dir = os.path.join(self.output_dirs[Utils.get_caller_name()], sampledata.name)
        assemblies = {"chromosome": "", "plasmid": ""}
        # spades is not defined in $PATH
        for assembly_type in assemblies:
            assembly_dir = os.path.join(stage_dir, assembly_type)
            cmd_append = ""
            if assembly_type == "plasmid":
                cmd_append = "--plasmid"
            cmd = f"""
            bash -c '
                cd {assembly_dir};
                export TOOL="$(find /usr/local/ -name "{_TOOL}*" -type f 2>/dev/null | grep '{_TOOL}.py$' | head -n 1)"; \
                echo "{_TOOL} is found at $TOOL";
                python3 "$TOOL" --version;
                python3 "$TOOL" \
                    --careful \
                    -o {assembly_dir} \
                    --threads {argValidator.threads} \
                    -1 {sampledata.reads[0]} \
                    -2 {sampledata.reads[1]} \
                    {cmd_append};
                chmod -R a+rw {assembly_dir};
            '
            """
            if not skip:
                self.clean_path(assembly_dir)
                log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
                Utils.append_log(log, _TOOL, sampledata.name)
            else:
                logging.info("Skip {} for {}".format(Utils.get_caller_name(), assembly_type))
            #  For most analyses, use scaffolds
            assembly_file = os.path.join(assembly_dir, "scaffolds.fasta")
            if not os.path.isfile(assembly_file):
                assembly_file = os.path.join(assembly_dir, "contigs.fasta")
                if not os.path.isfile(assembly_file):
                    logging.warning("The genome assemblies are missing for the sample: '{}'".format(sampledata.name))
            assemblies[assembly_type] = assembly_file
        sampledata.chromosome_assembly = assemblies["chromosome"]
        sampledata.plasmid_assembly = assemblies["plasmid"]

    def run_plasmid_merger(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "merge_chromosome_and_plasmid_assemblies"
        """
        # Sample launch:
        export IMG=ivasilyev/curated_projects:latest && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_dir = os.path.join(self.output_dirs[Utils.get_caller_name()], sampledata.name)
        genome_assembly = os.path.join(stage_dir, "{}_genome.fna".format(sampledata.name))
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 ./meta/scripts/merge_chromosome_and_plasmid_assemblies.py \
                --chromosome {sampledata.chromosome_assembly} \
                --plasmid {sampledata.plasmid_assembly} \
                --output {genome_assembly};
        '
        """
        if not skip:
            self.clean_path(stage_dir)
            log = Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip {}".format(Utils.get_caller_name()))
        if Utils.is_file_valid(genome_assembly):
            sampledata.genome_assembly = genome_assembly

    @staticmethod
    def _parse_taxa_from_blast_result(s: str):
        taxa_part = s.split("| ")[-1]
        d = Utils.parse_taxa(taxa_part)
        if len(d.keys()) < 0:
            logging.warning(f"Cannot parse taxa string: '{s}' -> '{d}'")
        return d

    def run_blast(self, sampledata: SampleDataLine, skip: bool = False):
        _TOOL = "blast_nucleotide_sequence"
        stage_name = Utils.get_caller_name()
        tool_dir = self.output_dirs[stage_name]
        stage_dir = os.path.join(tool_dir, sampledata.name)
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 ./meta/scripts/blast_nucleotide_sequence.py \
                --input {sampledata.genome_assembly} \
                --chromosomes_only \
                --results {argValidator.blast_reference_number} \
                --sequence_dir {self.blast_reference_dir} \
                --output {stage_dir};
        '
        """
        if not skip:
            self.clean_path(stage_dir)
            log = Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)
            Utils.append_log(log, _TOOL, sampledata.name)
        else:
            logging.info("Skip {}".format(Utils.get_caller_name()))
        blast_result_jsons = [
            i for i in Utils.scan_whole_dir(stage_dir)
            if os.path.basename(i).startswith(sampledata.name)
            and os.path.basename(i).endswith("_blast_results.json")
        ]
        if len(blast_result_jsons) == 0:
            logging.warning("No BLAST JSON result!")
            return
        sampledata.blast_result_json = blast_result_jsons[0]
        blast_result_dict = Utils.load_dict(sampledata.blast_result_json)
        blast_top_result = list(blast_result_dict.keys())[0]
        logging.info(f"The best matching organism is '{blast_top_result}'")

        taxa_dict = {i: "" for i in ["genus", "species", "strain"]}
        for blast_result in blast_result_dict.keys():
            taxa_dict.update(self._parse_taxa_from_blast_result(blast_result))
            if taxa_dict["species"] != "sp.":
                break
            logging.info(f"Unspecified species, trying another: '{blast_result}'")
        logging.info("The matching organism is {} {}".format(taxa_dict["genus"], taxa_dict["species"]))
        sampledata.set_taxa(**taxa_dict)

        sampledata.blast_result_table = Utils.locate_file_by_tail(stage_dir, "combined_blast_results.tsv")

        blast_report_json = Utils.locate_file_by_tail(stage_dir, "report.json")
        if len(blast_report_json) == 0:
            if not skip:
                logging.warning("No BLAST JSON report!")
            return
        blast_report_dict = Utils.load_dict(blast_report_json)
        genbank_files = blast_report_dict.get("genbank_files")
        if genbank_files is None or len(genbank_files) == 0:
            if not skip:
                logging.warning(f"Invalid BLAST JSON report: '{blast_report_json}'")
            return
        sampledata.closest_reference_genbank = genbank_files[0]
        if len(sampledata.closest_reference_genbank) == 0:
            logging.warning("No valid closest reference GenBank file is found!")
        self.update_state({
            sampledata.name: {
                stage_name: dict(
                    top_result_name=blast_top_result,
                    top_result_genbank=sampledata.closest_reference_genbank,
                    taxa=taxa_dict,
                )
            }
        })

    @staticmethod
    def _parse_quast_result(directory: str):
        """
        :param directory: str

        QUAST result directory

        :return: dict

        There are usually multiple report files with similar names:
        './report.tsv', './transposed_report.tsv', './reads_stats/reads_report.tsv'
        'report.tsv' contains two columns: key column and value column for the each assembly file.
        The key column:
        
        Assembly
        # contigs (>= 0 bp)
        # contigs (>= 1000 bp)
        # contigs (>= 5000 bp)
        # contigs (>= 10000 bp)
        # contigs (>= 25000 bp)
        # contigs (>= 50000 bp)
        Total length (>= 0 bp)
        Total length (>= 1000 bp)
        Total length (>= 5000 bp)
        Total length (>= 10000 bp)
        Total length (>= 25000 bp)
        Total length (>= 50000 bp)
        # contigs
        Largest contig
        Total length
        GC (%)
        N50
        N75
        L50
        L75
        # total reads
        # left
        # right
        Mapped (%)
        Properly paired (%)
        Avg. coverage depth
        Coverage >= 1x (%)
        # N's per 100 kbp
        """
        report_tables = sorted(Utils.locate_file_by_tail(directory, "report.tsv", multiple=True),
                               key=len,
                               reverse=False)
        out = dict()
        if len(report_tables) == 0:
            return out
        arr = Utils.load_2d_array(report_tables[0])
        for i in arr:
            if len(i) > 1:
                if i[1].isnumeric():
                    pass
                out[i[0]] = i[1]
        return out

    def run_quast(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "quast"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/quast:5.0.2--py36pl5262h30a8e3e_4 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_name = Utils.get_caller_name()
        stage_dir = os.path.join(self.output_dirs[stage_name], sampledata.name)

        genome_assembly_extension = Utils.get_file_extension(sampledata.genome_assembly)
        genome_assembly_symlink = os.path.join(stage_dir, f"{os.path.basename(sampledata.name)}{genome_assembly_extension}")

        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd "{stage_dir}";
            {_TOOL} \
                --gene-finding \
                --no-gzip \
                --output-dir "{stage_dir}" \
                --pe1 "{sampledata.raw_reads[0]}" \
                --pe2 "{sampledata.raw_reads[1]}" \
                --plots-format "png" \
                --features "{sampledata.closest_reference_genbank}" \
                --threads {argValidator.threads} \
                "{genome_assembly_symlink}"
            chmod -R a+rw "{stage_dir}";
        '
        """

        if skip:
            logging.info("Skip {}".format(stage_name))
        else:
            if len(sampledata.closest_reference_genbank) == 0:
                logging.warning("No valid closest reference GenBank file is found, aborting")
                return
            self.clean_path(stage_dir)
            os.symlink(sampledata.genome_assembly, genome_assembly_symlink)
            log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
            Utils.append_log(log, _TOOL, sampledata.name)
        quast_results = self._parse_quast_result(stage_dir)
        quast_result_number = len(quast_results.keys())
        if quast_result_number > 0:
            self.update_state({
                sampledata.name: {
                    stage_name: quast_results
                }
            })
        else:
            if not skip:
                logging.warning("No QUAST results!")

    @staticmethod
    def _locate_annotated_genome(directory: str):
        d = {i: "" for i in ["gbk", "gff"]}
        for extension in d.keys():
            files = [i for i in Utils.scan_whole_dir(directory) if i.endswith(".{}".format(extension))]
            if len(files) > 0:
                d[extension] = files[0]
            else:
                logging.warning("No annotated genome of type '{}'".format(extension))
        return d

    def run_prokka(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "prokka"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/prokka:1.12--pl526_0 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        stage_dir = os.path.join(self.output_dirs[Utils.get_caller_name()], sampledata.name)
        sampledata.genbank = os.path.join(stage_dir, "{}.gbf".format(sampledata.name))
        sampledata.faa = os.path.join(stage_dir, "{}.faa".format(sampledata.name))
        if skip or not sampledata.is_valid:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        taxa_append = ""
        for taxon_name, taxon_value in zip(
            ["genus", "species", "strain"],
            [sampledata.taxa_genus, sampledata.taxa_species, sampledata.taxa_strain]
        ):
            if taxon_value is not None and len(taxon_value) > 0:
                taxa_append = f"{taxa_append} --{taxon_name} {taxon_value}"
        if len(taxa_append) == 0:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        self.clean_path(stage_dir)

        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd {stage_dir};
            {_TOOL} \
                --centre UoN \
                --compliant \
                --cpu {argValidator.threads} \
                --force \
                --locustag {sampledata.name} \
                --outdir {stage_dir} \
                --prefix {sampledata.name} \
                --rfam \
                {taxa_append} \
                {sampledata.genome_assembly};
            chmod -R a+rw {stage_dir}
        '
        """
        log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
        Utils.append_log(log, _TOOL, sampledata.name)
        sampledata.genomes.update(self._locate_annotated_genome(stage_dir))

        gff_file = Utils.locate_file_by_tail(stage_dir, ".gff")
        if len(gff_file) == 0:
            return
        os.makedirs(self.roary_reference_dir, exist_ok=True)
        copy2(gff_file, os.path.join(self.roary_reference_dir, os.path.basename(gff_file)))

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

    def _run_getmlst(self, taxa: str, out_dir: str, out_file: str, sample_name: str):
        # One per sample, full cleaning is NOT required
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/srst2:0.2.0--py27_2 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
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
        _GETMLST_ATTEMPTS = 5
        os.makedirs(out_dir, exist_ok=True)
        for getmlst_attempt in range(1, _GETMLST_ATTEMPTS + 1):
            getmlst_cmd = """
            bash -c \
                'cd {o};
                 getmlst.py --species "{t}";
                 chmod -R a+rw {o}'
            """.format(t=taxa, o=out_dir)
            getmlst_log = self.run_quay_image("srst2", cmd=getmlst_cmd)
            Utils.append_log(getmlst_log, "getmlst", sample_name)
            srst2_cmd = getmlst_log.split("Suggested srst2 command for use with this MLST database:")[-1].strip()
            if not srst2_cmd.startswith("srst2"):
                logging.warning("`getmlst.py` did not finish correctly for attempt {} of {}".format(
                    getmlst_attempt, _GETMLST_ATTEMPTS))
            else:
                logging.info("Got the SRST2 output: '{}'".format(srst2_cmd))
                out = {j[0]: " ".join(j[1:]) for j in [i.strip().split(" ") for i in srst2_cmd.split("--")[1:]]}
                out.update({i: os.path.join(out_dir, out[i]) for i in ["mlst_db", "mlst_definitions"]})
                Utils.dump_dict(out, out_file)
                return out
            if getmlst_attempt == _GETMLST_ATTEMPTS:
                logging.warning("Exceeded attempts number for `getmlst.py` to finish processing correctly")

    @staticmethod
    def _parse_srst2_result_log(srst2_out_mask: str):
        # Extract information from the internal SRST2 log
        _PHRASE = " MLST output printed to "
        log_file = "{}.log".format(srst2_out_mask)
        if Utils.is_file_valid(log_file):
            log_lines = [j for j in Utils.remove_empty_values(
                [i.strip() for i in Utils.load_list(log_file)]) if _PHRASE in j]
            if len(log_lines) > 0:
                srst2_result_file = log_lines[0].split(_PHRASE)[-1].strip()
                return srst2_result_file

    # MLST typing
    def run_srst2(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample, full cleaning is NOT required
        _TOOL = "srst2"
        _SRST2_ATTEMPTS = 5
        _ERROR_PHRASES = [
            "Encountered internal Bowtie 2 exception",  "[main_samview] truncated file."
        ]
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/srst2:0.2.0--py27_2 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        reference_dir = os.path.join(self.srst2_reference_dir, sampledata.prefix)
        stage_dir = os.path.join(tool_dir, sampledata.name)
        if skip or sampledata.is_taxa_valid is None:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        self.clean_path(stage_dir)
        os.makedirs(stage_dir, exist_ok=True)
        genus, species = sampledata.taxa_genus, sampledata.taxa_species
        taxa = " ".join([genus, species])
        getmlst_state_file = os.path.join(reference_dir, "getmlst_state_{}.json".format(taxa))
        if Utils.is_file_valid(getmlst_state_file):
            getmlst_state = Utils.load_dict(getmlst_state_file)
        else:
            getmlst_state = self._run_getmlst(taxa=taxa,
                                              out_dir=reference_dir,
                                              out_file=getmlst_state_file,
                                              sample_name=sampledata.name)
        """
        SRST2 requires the strict file system structure pattern:

        <work_directory>/
            ├── <sample_name>_1.<reads_extension>
            ├── <sample_name>_2.<reads_extension>
            ├── <taxonomy_prefix>.log
            ├── <taxonomy_prefix>__<sample_name>.<Genera>_<species>.pileup
            ├── <taxonomy_prefix>__<sample_name>.<Genera>_<species>.sorted.bam
            └── <taxonomy_prefix>__mlst__<Genera>_<strain>__results.txt
         <reference_directory>/
            ├── alleles_fasta
            ├── getmlst_state_<Genera>_<species>_<strain>.json
            ├── profiles_csv
            ├── <Genera>_<strain>.fasta
            ├── <Genera>_<strain>.fasta.<bowtie_2_index_number>.bt2
            ├── <Genera>_<strain>.fasta.<bowtie_2_index_number>.bt2
            ├── <Genera>_<strain>.fasta.fai
            └── <Genera>_<strain>.fasta.rev.<bowtie_2_rev_index_number>.bt2

        More info: https://github.com/katholt/srst2#input-read-formats-and-options
        """

        input_reads = [os.path.join(stage_dir, f"{sampledata.name}_{idx + 1}{sampledata.extension}")
                       for idx, i in enumerate(sampledata.reads)]
        out_mask = os.path.join(stage_dir, sampledata.prefix)

        for file, link in zip(sampledata.reads, input_reads):
            os.symlink(file, link)

        cmd = f"""
        bash -c '
            {_TOOL} --version;
            cd {stage_dir};
            {_TOOL} --log \
               --input_pe {" ".join(input_reads)} \
               --mlst_db {getmlst_state["mlst_db"]} \
               --mlst_definitions {getmlst_state["mlst_definitions"]} \
               --mlst_delimiter {getmlst_state["mlst_delimiter"]} \
               --output {out_mask} \
               --threads {argValidator.threads};
            chmod -R a+rw {stage_dir}
        '
        """
        # Deliberately set the tag with fully supported environment
        # https://quay.io/repository/biocontainers/srst2?tab=tags
        log = self.run_quay_image(
            _TOOL, img_tag="0.2.0--py27_2", cmd=cmd, attempts=_SRST2_ATTEMPTS,
            bad_phrases=_ERROR_PHRASES, sample_name=sampledata.name
        )
        Utils.append_log(log, _TOOL, sampledata.name)

        sampledata.srst2_result_table = self._parse_srst2_result_log(out_mask)
        if not os.path.isfile(sampledata.srst2_result_table):
            logging.warning("Not found the SRST2 processing result file: '{}', trying to locate it".format(
                sampledata.srst2_result_table))
            sampledata.srst2_result_table = Utils.locate_file_by_tail(stage_dir, "__results.txt")

    def merge_srst2_results(self, sampledata_array: SampleDataArray, skip: bool = False):
        _TOOL = "concatenate_tables"
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        if skip:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        sampledata_array.srst2_merged_table = os.path.join(tool_dir, "srst2_merged_results.tsv")
        table_list = sampledata_array.srst2_result_tables
        if len(table_list) == 0:
            logging.warning("No SRST@ result tables!")
            return
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 ./meta/scripts/{_TOOL}.py \
                --input {Utils.render_file_list(table_list)} \
                --axis 0 \
                --output {sampledata_array.srst2_merged_table};
        '
        """
        self.clean_path(tool_dir)
        log = Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)
        Utils.append_log(log, _TOOL)

    def _download_card_reference(self):
        # The CARD reference updates relatively frequent, so it's better to fetch a fresh copy per launch
        self.clean_path(self.card_reference_dir)
        cmd = f"""
        bash -c '
            cd {self.card_reference_dir};
            for i in "https://card.mcmaster.ca/latest/data" "https://card.mcmaster.ca/latest/ontology" "https://card.mcmaster.ca/latest/variants";
                do \
                    curl -fsSL "$i" | tar jxf -;
                done;
        '
        """
        logging.info("Download the CARD reference")
        return Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)

    def run_rgi(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "rgi"
        """
        # Sample launch:
        export IMG=quay.io/biocontainers/rgi:5.1.1--py_0 && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        stage_dir = os.path.join(tool_dir, sampledata.name)
        if skip or not sampledata.is_valid:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return

        self.card_reference_json = argValidator.card_json
        if not Utils.is_file_valid(self.card_reference_json, report=True):
            os.makedirs(self.card_reference_dir, exist_ok=True)
            for _try in range(0, PULL_RETRIES):
                self.card_reference_json = Utils.locate_file_by_tail(self.card_reference_dir, "card.json")
                if len(self.card_reference_json) > 0:
                    break
                logging.info(f"Downloading CARD reference data for attempt {_try} of {PULL_RETRIES}")
                log = self._download_card_reference()
                Utils.append_log(log, _TOOL, sampledata.name)

        self.clean_path(stage_dir)
        out_mask = os.path.join(stage_dir, sampledata.name)

        # There is no '-v/--version' CLI argument
        # Outputs here are masks only
        cmd = f"""
        bash -c '
            echo "$({_TOOL} --help | grep '^Resistance')";
            cd "{stage_dir}";
            {_TOOL} load --card_json "{self.card_reference_json}";
            {_TOOL} main \
                --clean \
                --input_sequence "{sampledata.genome_assembly}" \
                --input_type contig \
                --num_threads {argValidator.threads} \
                --output_file "{out_mask}";
            for category in "drug_class" "resistance_mechanism" "gene_family";
                do \
                    {_TOOL} heatmap \
                        --input "{stage_dir}" \
                        --category "$category" \
                        --output "{out_mask}_heatmap_by_$category";
                done;
            chmod -R a+rw "{stage_dir}";
        '
        """
        log = self.run_quay_image(_TOOL, cmd=cmd, sample_name=sampledata.name)
        Utils.append_log(log, _TOOL, sampledata.name)
        """
        Heatmap info:
        
        Yellow represents a perfect hit, 
        teal represents a strict hit, 
        purple represents no hit. 
        Genes with asterisks (*) appear multiple times.
        """

    @staticmethod
    def _convert_and_index_reference_by_bwa(input_genbank_file: str, output_fna_file: str):
        # One per sample
        """
        # Sample launch:
        export IMG="ivasilyev/mgefinder:latest" && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}"
        """
        output_dir = os.path.dirname(output_fna_file)
        cmd = f"""
        bash -c '
            cd {output_dir} && \
            python3 "$CONVERTER" \
                --input "{input_genbank_file}" \
                --input_format genbank \
                --output "{output_fna_file}" \
                --output_format fasta;
            bwa index \
                -p "{Utils.get_file_name_mask(output_fna_file)}" \
                "{output_fna_file}";
            chmod -R a+rw {output_dir};
        '
        """
        return Utils.run_image(img_name="ivasilyev/mgefinder:latest", container_cmd=cmd)

    def run_mgefinder(self, sampledata: SampleDataLine, skip: bool = False):
        # One per sample
        _TOOL = "mgefinder"
        """
        # Sample launch:
        export IMG="ivasilyev/mgefinder:latest" && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}"
        """
        stage_name = Utils.get_caller_name()
        tool_dir = self.output_dirs[stage_name]
        stage_dir = os.path.join(tool_dir, sampledata.name)

        if skip or not sampledata.is_valid:
            logging.info("Skip {}".format(stage_name))
            return

        os.makedirs(self.mgefinder_reference_dir, exist_ok=True)

        """
        MGEfinder requires the very strict file system structure pattern. 
        File extensions are affected as well: 
        
        <work_directory>/
            ├── 00.assembly/
            │   ├── <sample_name>.fna
            │   ...
            ├── 00.bam/
            │   ├── <sample_name>.<reference_name>.bam
            │   ├── <sample_name>.<reference_name>.bam.bai
            │   ...
            └── 00.genome/
                └── <reference_name>.fna
        
        More info: https://github.com/bhattlab/MGEfinder/wiki
        """

        assembly_dir = os.path.join(stage_dir, "00.assembly")
        bam_dir = os.path.join(stage_dir, "00.bam")
        genome_dir = os.path.join(stage_dir, "00.genome")
        raw_dir = os.path.join(stage_dir, "00.raw")
        _ = [self.clean_path(i) for i in (stage_dir, assembly_dir, bam_dir, genome_dir, raw_dir)]

        reads_file_1 = os.path.join(raw_dir, f"{sampledata.name}.1{sampledata.extension}")
        reads_file_2 = os.path.join(raw_dir, f"{sampledata.name}.2{sampledata.extension}")

        assembly_name = os.path.join(assembly_dir, f"{sampledata.name}.fna")

        if len(sampledata.closest_reference_genbank) == 0:
            logging.warning("No valid closest reference GenBank file is found, aborting")
            return
        reference_name = Utils.filename_only(sampledata.closest_reference_genbank)
        genome_file = f"{os.path.join(genome_dir, reference_name)}.fna"
        genome_mask = Utils.get_file_name_mask(genome_file)

        if not Utils.is_file_valid(genome_file):
            logging.info("Create BWA index")
            log = self._convert_and_index_reference_by_bwa(
                sampledata.closest_reference_genbank, genome_file
            )
            Utils.append_log(log, "bwa index")

        alignment_mask = os.path.join(bam_dir, f"{sampledata.name}.{reference_name}")

        for file, link in (
            [sampledata.reads[0], reads_file_1],
            [sampledata.reads[1], reads_file_2],
            [sampledata.genome_assembly, assembly_name]
        ):
            os.symlink(file, link)

        cmd = f"""
        bash -c '
            mgefinder --version;
            echo "Align sequences";
            bwa mem \
                -t {argValidator.threads} \
                -v 0 \
                "{genome_mask}" \
                "{reads_file_1}" \
                "{reads_file_2}" \
                > "{alignment_mask}.sam";
            echo "Convert SAM to BAM";
            mgefinder formatbam \
                "{alignment_mask}.sam" \
                "{alignment_mask}.bam";
        '
        """
        log = Utils.run_image(img_name="ivasilyev/mgefinder:latest", container_cmd=cmd)
        Utils.append_log(log, _TOOL)

        cmd = f"""
        bash -c '
            cd "{stage_dir}";
            echo "Start MGEfinder";
            mgefinder workflow denovo \
                --cores {argValidator.threads} \
                "{stage_dir}";
            chmod -R a+rw "{stage_dir}";
        '
        """
        log = Utils.run_image(img_name="ivasilyev/mgefinder:latest", container_cmd=cmd)
        Utils.append_log(log, _TOOL)

    # Group processing pipeline steps

    def merge_blast_results(self, sampledata_array: SampleDataArray, skip: bool = False):
        _TOOL = "concatenate_tables"
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        if skip:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        sampledata_array.blast_merged_table = os.path.join(tool_dir, "blast_merged_results.tsv")
        table_list = sampledata_array.blast_result_tables
        if len(table_list) == 0:
            logging.warning("No BLAST result tables!")
            return
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 ./meta/scripts/{_TOOL}.py \
                --input {Utils.render_file_list(table_list)} \
                --axis 0 \
                --output {sampledata_array.blast_merged_table};
        '
        """
        self.clean_path(tool_dir)
        log = Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)
        Utils.append_log(log, _TOOL)

    @staticmethod
    def _convert_genbank_to_gff3(gbff_dir, gff_dir):
        # One per all samples, force multithread
        _TOOL = "bp_genbank2gff3"
        """
        # Sample launch:
        export IMG=bioperl/bioperl:latest && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash

        # Tool executable lookup:
        export TOOL="$(find /usr/local/ -name "bp_genbank2gff3" -type f 2>/dev/null | grep 'bp_genbank2gff3$' | head -n 1)"
        echo "${TOOL}"
        """
        cmd_1 = f"""
        export TOOL="$(find /usr/local/ -name "{_TOOL}" -type f 2>/dev/null | grep '{_TOOL}' | head -n 1)";
        export SOURCE_DIR="{gbff_dir}/";
        export TARGET_DIR="{gff_dir}/";
        """
        cmd_2 = """
        cd "${TARGET_DIR}";
        find "${SOURCE_DIR}" \
            -name '*.gbk' \
            -type f \
            -print0 \
        | xargs \
            -0 \
            --max-procs "$(nproc)" \
            -I "{}" \
            bash -c '
                TARGET_FILE="${TARGET_DIR}$(basename "{}").gff";
                if [ ! -s "${TARGET_FILE}" ]; 
                    then
                        perl "${TOOL}" \
                            "{}" \
                            --outdir "${TARGET_DIR}";
                fi;
            ';
        chmod -R a+rw "${TARGET_DIR}";
        """
        # Mock the script, since the `bash -c 'bash -c '...''` hangs
        exe = os.path.join(gbff_dir, f"{_TOOL}.sh")
        if not Utils.is_file_valid(exe):
            Utils.dump_string(cmd_2, exe)
        # Roary fails if found a non-GFF file in its input
        cmd = f"""
        bash -c '{cmd_1}
        bash {exe};
        '
        """
        os.makedirs(gff_dir, exist_ok=True)
        return Utils.run_image(img_name="bioperl/bioperl:latest", container_cmd=cmd)

    @staticmethod
    def _process_newick(newick_file: str, blast_result_table: str, out_file: str):
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 ./meta/scripts/replace_ids_with_strains_from_blast_result_table.py \
                --input {newick_file} \
                --table {blast_result_table} \
                --output {out_file};
        '
        """
        return Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)

    def run_roary(self, sampledata_array: SampleDataArray, skip: bool = False):
        # One per all samples
        _TOOL = "roary"
        """
        # Sample launch:
        export IMG=sangerpathogens/roary:latest && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        stage_dir = os.path.join(tool_dir, "out")
        # Stage dir here must be created by the program only
        if skip:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        self.clean_path(tool_dir)
        os.makedirs(self.roary_reference_dir, exist_ok=True)
        log = self._convert_genbank_to_gff3(self.blast_reference_dir, self.roary_reference_dir)
        Utils.append_log(log, _TOOL)

        # There is no version number accessible via special or even help CLI arguments
        cmd = f"""bash -c '
            cd {tool_dir};
            {_TOOL} \
                -f "{stage_dir}" \
                -e \
                --mafft \
                -p {argValidator.threads} \
                "{self.roary_reference_dir}"/*;
            chmod -R a+rw "{tool_dir}";
        '
        """
        log = Utils.run_image(img_name="sangerpathogens/roary:latest", container_cmd=cmd)
        Utils.append_log(log, _TOOL)

        sampledata_array.roary_edited_newick = os.path.join(tool_dir, "roary_edited_results.newick")
        newick_file = Utils.locate_file_by_tail(stage_dir, ".newick")
        if len(newick_file) > 0:
            logging.info(f"Found Newick file: '{newick_file}'")
            log = self._process_newick(
                newick_file,
                sampledata_array.blast_merged_table,
                sampledata_array.roary_edited_newick
            )
            Utils.append_log(log, _TOOL)
        else:
            logging.warning("No Newick file found!")

    @staticmethod
    def _run_nbee(sampledata_file: str, refdata_file: str, out_dir: str):
        # One per all samples
        _TOOL = "nBee"
        """
        # Sample launch:
        export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        cmd = f"""
        bash -c '
            python3 "$HOME/scripts/{_TOOL}.py" \
                --input "{sampledata_file}" \
                --refdata "{refdata_file}" \
                --output "{out_dir}"
        '
        """
        log = Utils.run_image(img_name="ivasilyev/bwt_filtering_pipeline_worker:latest",
                              container_cmd=cmd)
        if len(log) > 0:
            Utils.append_log(log, _TOOL)

    @staticmethod
    def _annotate_and_filter_nbee_coverage(coverage_table: str, annotation_table: str,
                                           annotated_table: str, filtered_table):
        # One per sample
        _TOOL = "concatenate_tables+filter_only_not_null_values"
        """
        # Sample launch:
        export IMG=ivasilyev/curated_projects:latest && \
        docker pull "${IMG}" && \
        docker run --rm --net=host -it "${IMG}" bash
        """
        cmd = f"""
        bash -c '
            git pull --quiet && \
            python3 "$HOME/scripts/curated_projects/meta/scripts/concatenate_tables.py" \
                --axis 1 \
                --index "reference_id" \
                --input "{annotation_table}" "{coverage_table}" \
                --join inner \
                --output "{annotated_table}" && \
            python3 "$HOME/scripts/curated_projects/meta/scripts/filter_only_not_null_values.py" \
                --input "{annotated_table}" \
                --filter "id_mapped_relative_abundance" \
                --output "{filtered_table}"
        '
        """
        log = Utils.run_image(img_name="ivasilyev/curated_projects:latest", container_cmd=cmd)
        if len(log) > 0:
            Utils.append_log(log, _TOOL)

    def _run_and_annotate_nbee(self, sampledata_file: str, refdata_file: str, tool_dir: str):
        try:
            refdata_dict = Utils.load_dict(refdata_file)
            refdata_first_dict = refdata_dict[list(refdata_dict.keys())[0]]
            annotation_file = refdata_first_dict["annotation"]
            refdata_alias = refdata_first_dict["alias"]
            out_dir = os.path.join(tool_dir, refdata_alias)
        except json.JSONDecodeError:
            logging.warning(f"Invalid RefData format for '{refdata_file}'")
            return

        # Mapping
        self._run_nbee(sampledata_file=sampledata_file, refdata_file=refdata_file, out_dir=out_dir)

        # Annotation
        coverage_files = Utils.locate_file_by_tail(out_dir, "_coverage.tsv", multiple=True)
        annotation_dir = os.path.join(out_dir, "annotated_coverages")
        os.makedirs(annotation_dir, exist_ok=True)
        for raw_coverage_file in coverage_files:
            coverage_file_extension = Utils.get_file_extension(raw_coverage_file)
            annotated_coverage_file = os.path.join(
                annotation_dir,
                "{}_annotated{}".format(
                    os.path.basename(os.path.splitext(raw_coverage_file)[0]),
                    coverage_file_extension
                )
            )
            filtered_coverage_file = "{}_filtered{}".format(
                os.path.splitext(annotated_coverage_file)[0],
                coverage_file_extension
            )
            self._annotate_and_filter_nbee_coverage(
                coverage_table=raw_coverage_file,
                annotation_table=annotation_file,
                annotated_table=annotated_coverage_file,
                filtered_table=filtered_coverage_file,
            )

    def run_nbee_with_annotation(self, sampledata_array: SampleDataArray, skip: bool = False):
        # One per all samples per all refdata
        _TOOL = "nBee"
        stage_name = Utils.get_caller_name()
        tool_dir = self.output_dirs[stage_name]
        self.clean_path(tool_dir)
        if skip:
            logging.info("Skip {}".format(stage_name))
            return

        sampledata_file = os.path.join(tool_dir, "sampledata.tsv")
        Utils.dump_2d_array(sampledata_array.to_2d_array(), sampledata_file)

        refdata_files = [i for i in argValidator.refdata_files if Utils.is_file_valid(i)]
        if len(refdata_files) == 0:
            logging.warning("No valid RefData files found, the alignment aborted!")
            return

        for refdata_file in refdata_files:
            self._run_and_annotate_nbee(sampledata_file=sampledata_file,
                                        refdata_file=refdata_file,
                                        tool_dir=tool_dir)

    # Orthologs-based phylogenetic tree construction
    def run_orthomcl(self, sampledata_array: SampleDataArray, skip: bool = False):
        return
        # One per all samples
        _TOOL = "orthomcl"
        tool_dir = self.output_dirs[Utils.get_caller_name()]
        if skip:
            logging.info("Skip {}".format(Utils.get_caller_name()))
            return
        self.clean_path(tool_dir)
        sampledata_file = os.path.join(tool_dir, "orthomcl_input.sampledata")
        logging.info("Save OrthoMCL sample data: '{}'".format(sampledata_file))
        Utils.dump_2d_array([[i.name, i.genbank] for i in sampledata_array.lines], sampledata_file)
        log = Utils.run_image(img_name="ivasilyev/orthomcl-mysql:latest",
                              container_cmd="""
                              bash -c \
                                'service mysql restart;
                                 python3 /opt/my_tools/pipeline_handler.py -i {s} -o {o};
                                 chmod -R a+rw {o}'
                              """.format(s=sampledata_file, o=tool_dir))
        Utils.append_log(log, _TOOL, "all")

    def dump_state(self):
        Utils.dump_dict(self._state, os.path.join(self.output_dir_root, "state.json"))

    def update_state(self, d: dict):
        Utils.merge_dicts(source=d, target=self._state)
        self.dump_state()

    def handle(self, sampledata_array: SampleDataArray):
        if not self.valid:
            return
        for idx, func in enumerate(self.sample_methods + self.group_methods):
            try:
                logging.info("Starting the pipeline step {} of {} ({} in total)".format(
                    idx + 1, len(self.sample_methods + self.group_methods), len(argValidator.stages_to_do)))
                # Per-sample processing
                if func in self.sample_methods:
                    queue = [(func, i, idx not in argValidator.stages_to_do) for i in sampledata_array.lines.values()]
                    _ = Utils.single_core_queue(Utils.wrap_func, queue)
                # Per-group functions
                elif func in self.group_methods:
                    _ = func(sampledata_array, skip=idx not in argValidator.stages_to_do)
            except PermissionError:
                logging.critical("Cannot process the step {}, please run the command 'sudo chmod -R a+rw {}'".format(
                    idx, self.output_dir_root))
        self.dump_state()


class Utils:
    # File system based methods
    @staticmethod
    def is_file_valid(file: str, report: bool = False):
        if len(file) == 0:  # Do not log this
            return False
        if not os.path.exists(file):
            if report:
                logging.warning("Not found: '{}'".format(file))
            return False
        if not os.path.isfile(file):
            if report:
                logging.warning("Not a file: '{}'".format(file))
            return False
        if os.path.getsize(file) == 0:
            if report:
                logging.warning("Empty file: '{}'".format(file))
            return False
        return True

    @staticmethod
    def scan_whole_dir(dir_name: str):
        out = []
        for root, dirs, files in os.walk(dir_name):
            for file in files:
                out.append(os.path.join(root, file))
        return sorted(out)

    @staticmethod
    def get_file_name_mask(path: str):
        return os.path.splitext(os.path.normpath(path))[0]

    @staticmethod
    def filename_only(path: str):
        return os.path.splitext(os.path.basename(os.path.normpath(path)))[0]

    @staticmethod
    def get_file_extension(file: str, deep: int = 1):
        split = os.path.basename(file.strip()).split(".")[::-1]
        out = []
        for sub in split:
            if 5 >= len(sub) > 1:
                out.append(str(sub))
            else:
                break
        return ".{}".format(".".join(out[:deep][::-1]))

    @staticmethod
    def get_reads_file_extension(file: str):
        out = Utils.get_file_extension(file)
        _EXTENSIONS = {
            "fastq.gz": ("fastq.gz", "fq.gz"),
            "fastq": ("fastq", "fq")
        }
        for replacement, keywords in _EXTENSIONS.items():
            if any(file.endswith(f".{i}") for i in keywords):
                out = f".{replacement}"
        return out

    @staticmethod
    def locate_file_by_tail(dir_name: str, tail: str, multiple: bool = False):
        files = [i for i in Utils.scan_whole_dir(dir_name) if i.endswith(tail)]
        if len(files) == 0:
            return ""
        if multiple:
            return files
        return files[0]

    # Date / time methods

    @staticmethod
    def get_time():
        now = datetime.now()
        output_list = []
        for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
            time_unit = str(time_unit)
            if len(time_unit) < 2:
                time_unit = '0' + time_unit
            output_list.append(time_unit)
        return '-'.join(output_list)

    @staticmethod
    def count_elapsed_seconds(t):
        return f"{perf_counter() - t :.3f} s."

    # I/O methods

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
    def load_dict(file: str):
        return json.loads(Utils.load_string(file))

    @staticmethod
    def dump_dict(d: dict, file: str, **kwargs):
        _kwargs = dict(indent=4, sort_keys=False)
        if len(kwargs.keys()) > 0:
            _kwargs.update(kwargs)
        return Utils.dump_string(json.dumps(d, **_kwargs), file)

    @staticmethod
    def load_list(file: str):
        return Utils.split_lines(Utils.load_string(file))

    @staticmethod
    def dump_list(lst: list, file: str):
        Utils.dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)

    @staticmethod
    def list_to_spaced_string(x: list):
        return " ".join([str(i) for i in x])

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utils.dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)

    @staticmethod
    def append_log(msg: str, tool_name: str, sample_name: str = "all"):
        logging.debug(msg)
        file = os.path.join(argValidator.log_dir, "{}_{}.log".format(tool_name, sample_name))
        with open(file, mode="a", encoding="utf-8") as f:
            f.write(msg + "\n")
            f.close()

    @staticmethod
    def unzip_archive(archive: str, extract_path: str, remove=True):
        import shutil
        os.makedirs(extract_path, exist_ok=True)
        shutil.unpack_archive(archive, extract_path)
        print("Extracting completed: '{}'".format(archive))
        if remove:
            os.remove(archive)
            print("Removed: '{}'".format(archive))

    # Primitive processing methods

    @staticmethod
    def safe_findall(pattern, string, idx: int = 0, verbose: bool = False):
        try:
            return re.findall(pattern, string)[idx]
        except IndexError:
            if verbose:
                print(f"Warning! Can't find the regex pattern '{pattern}' within the string: '{string}'")
            return ""

    @staticmethod
    def split_lines(string: str):
        out = [i.strip() for i in re.sub("[\r\n]+", "\n", string).split("\n")]
        return Utils.remove_empty_values(out)

    @staticmethod
    def string_to_2d_array(string: str):
        out = [[j.strip() for j in i.split("\t")] for i in Utils.split_lines(string)]
        return Utils.remove_empty_values(out)

    @staticmethod
    def load_2d_array(file: str):
        return Utils.string_to_2d_array(Utils.load_string(file))

    @staticmethod
    def flatten_2d_array(array: list):
        return [j for i in array for j in i]

    @staticmethod
    def remove_empty_values(input_list):
        output_list = []
        if input_list is not None:
            for i in input_list:
                if i is not None:
                    try:
                        if len(i) > 0:
                            output_list.append(i)
                    except TypeError:
                        continue
        return output_list

    @staticmethod
    def parse_taxa(taxa):
        out = dict(genus="", species="", strain="")
        if isinstance(taxa, str):
            taxa_list = Utils.remove_empty_values(re.split("[. ]+", taxa.strip()))
            if len(taxa_list) == 0:
                return out
            out["genus"] = taxa_list[0]
            if len(taxa_list) > 1:
                if any(i.isdigit() for i in taxa_list[1]):  # Case 'Escherichia O157:H7'
                    out["strain"] = taxa_list[1]
                else:  # Case 'Escherichia coli'
                    sp = taxa_list[1]
                    if sp in "spp":  # Case 'Escherichia sp.'
                        out["species"] = "sp."
                    else:
                        out["species"] = sp
            if len(taxa_list) > 2:  # Case 'Escherichia coli O157:H7'
                out["strain"] = taxa_list[2]
        if isinstance(taxa, dict):
            out["genus"], out["species"], out["strain"] = [
                j if j is not None else ""
                for j in [
                    taxa.get(i) for i in out.keys()
                ]
            ]
        return out

    @staticmethod
    def render_file_list(x):
        return '"{}"'.format('", "'.join(Utils.remove_empty_values(x)))

    @staticmethod
    def merge_dicts(source: dict, target: dict):
        """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
        updating only top-level keys, dict_merge recurses down into dicts nested
        to an arbitrary depth, updating keys. The ``source`` is merged into
        ``target``.
        :param target: dict onto which the merge is executed
        :param source: target merged into target
        :return: None
        """
        for k, v in source.items():
            if k in target and isinstance(target[k], dict) and isinstance(source[k], Mapping):
                Utils.merge_dicts(source=source[k], target=target[k])
            else:
                target[k] = source[k]

    # Function handling methods

    @staticmethod
    def randomize_sleep(min_: int = 30, max_: int = 120):
        from time import sleep
        from random import randint
        sleep(randint(min_, max_))

    @staticmethod
    def get_caller_name():
        import inspect
        return str(inspect.stack()[1][3])

    @staticmethod
    def log_and_raise(msg):
        logging.critical(msg)
        raise ValueError(msg)

    @staticmethod
    def is_log_valid(log: str, bad_phrases: list):
        log_lines = Utils.remove_empty_values(log.replace("\r", "\n").split("\n"))
        bad_phrases = Utils.remove_empty_values(bad_phrases)
        if len(bad_phrases) == 0:
            return True
        return all([i not in j for i in bad_phrases for j in log_lines])

    @staticmethod
    def run_until_valid_output(
            cmd: str, bad_phrases: list, attempts: int = 5, ping_required: bool = False
    ):
        start = perf_counter()
        attempt = 0
        log = ""
        while attempt < attempts:
            attempt += 1
            log_cmd = re.sub("[\r\n ]+", " ", cmd.strip())
            logging.debug("Executing the command: `{}`".format(log_cmd))
            log = subprocess.getoutput(cmd)
            if not Utils.is_log_valid(log, bad_phrases):
                logging.warning(
                    f"An error phrase was found in log for attempt {attempt} of {attempts}."
                )
                sleep(5)
                if ping_required:
                    # Ping Google just to keep the node DNS working
                    _ = subprocess.getoutput("ping -c 10 google.com")
            else:
                break
            if attempt == attempts:
                logging.warning(
                    f"Exceeded attempts number to get execution output without failure messages. The command seems to be not finished correctly: `{log_cmd}`"
                )
        logging.debug(
            f"Completed after {attempt} attempt(s) and {Utils.count_elapsed_seconds(start)}"
        )
        return log

    @staticmethod
    def run_image(
            img_name: str, container_cmd: str, bad_phrases: list = (), attempts: int = 5
    ):
        _DOCKER_RUN_CMD = "docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it"
        _COMMON_PHRASES = ["Error response from daemon", ]
        logging.info("Using image: '{}'".format(img_name))
        bad_phrases = list(bad_phrases) + _COMMON_PHRASES
        docker_cmd = f'docker pull "{img_name}" && {_DOCKER_RUN_CMD} "{img_name}"'
        # cmd may contain curly braces, so str.format() is not usable
        out_cmd = docker_cmd.strip() + " " + container_cmd.strip()
        return Utils.run_until_valid_output(
            cmd=out_cmd, bad_phrases=bad_phrases, attempts=attempts, ping_required=True
        )

    @staticmethod
    def wrap_func(args: list):
        return args[0](*args[1:])

    # Queue processing methods

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

    # Web-based methods

    @staticmethod
    def download_file(url, out_file, force: bool = False):
        import subprocess
        from time import sleep
        _RETRIES = 5
        _SLEEP_SECONDS = 3
        _ERROR_REPORTS = ("transfer closed with",)
        url = url.strip()
        out_file = out_file.strip()
        if not force and Utils.is_file_valid(out_file):
            logging.debug("Already downloaded: '{}'".format(out_file))
            return
        out_dir = os.path.dirname(out_file)
        assert len(url) > 0 and len(out_dir) > 0
        os.makedirs(out_dir, exist_ok=True)
        for c in range(_RETRIES):
            log = subprocess.getoutput("curl -fsSL {} -o {}".format(url, out_file))
            print(log)
            if Utils.is_file_valid(out_file) and all(i not in log for i in _ERROR_REPORTS):
                logging.debug("Download finished: '{}'".format(out_file))
                return
            sleep(_SLEEP_SECONDS)
            logging.warning("Warning! Failed download: '{}'. Retries left: {}".format(url, _RETRIES - c - 1))
        logging.warning("Exceeded URL download limits: '{}'".format(url))

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


if __name__ == '__main__':
    mainStart = perf_counter()
    argValidator = ArgValidator()
    mainLogFile = os.path.join(argValidator.log_dir, "main.log")
    os.makedirs(argValidator.log_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[logging.FileHandler(mainLogFile), logging.StreamHandler()],
        format="asctime=%(asctime)s levelname=%(levelname)s process=%(process)d name=%(name)s "
               "funcName=%(funcName)s lineno=%(lineno)s message=\"\"\"%(message)s\"\"\""
    )
    argValidator.validate()
    sampleDataArray = SampleDataArray.load(argValidator.sampledata_file)
    handler = Handler(argValidator.output_dir)
    logging.info("The pipeline processing started")
    handler.handle(sampleDataArray)
    logging.info(f"The pipeline processing completed after {Utils.count_elapsed_seconds(mainStart)}")
