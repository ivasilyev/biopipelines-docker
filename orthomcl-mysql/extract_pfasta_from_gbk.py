#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class ArgParser:
    def __init__(self):
        import argparse
        _parser = argparse.ArgumentParser(description="Extract amino acid sequence from provided GenBank file.".strip(),
                                          epilog="Note: the input file must be compliant with actual GenBank standard, "
                                                 "otherwise it won't be parsed with BioPython")
        _parser.add_argument("-i", "--input", metavar="<input.gbk>", type=str, required=True, help="Input GenBank file")
        _parser.add_argument("-a", "--abbreviation", metavar="<str>", type=str, default="",
                             help="Organism species abbreviation containing 3 or 4 characters")
        _parser.add_argument("-s", "--sample_name", metavar="<str>", type=str, default="",
                             help="Sample name to add into prefix")
        _parser.add_argument("-o", "--output", metavar="<output.faa>", type=str, required=True, help="Output file name")
        self._namespace = _parser.parse_args()
        self.input_gbk = self._namespace.input
        self.abbreviation = self._namespace.abbreviation.capitalize()
        self.sample_name = self._namespace.sample_name
        self.out_pfasta = self._namespace.output


class Converter:
    def __init__(self):
        self.input_gbk = parser.input_gbk
        self.abbreviation = parser.abbreviation
        self.sample_name = parser.sample_name
        self.out_pfasta = parser.out_pfasta
        self._out_pfasta_records = []
        self._out_annotations = dict()

    def parse(self):
        seq_records = list(SeqIO.parse(self.input_gbk, "genbank"))
        _id = 0
        # LOCUS
        for seq_record in seq_records:
            # CDS
            for seq_feature in seq_record.features:
                qualifiers = seq_feature.qualifiers.copy()
                if seq_feature.type == "CDS" and "translation" in qualifiers.keys():
                    _id += 1
                    id_str = "PFASTA_ID_{}".format(_id)
                    contig = "contig_{}".format(Utils.safe_extract_int(seq_record.id))
                    locus_tag = "CDS_{}".format(Utils.safe_extract_int(Utils.safe_get(qualifiers, "locus_tag")))
                    gene = Utils.safe_get(qualifiers, "gene")
                    product = Utils.safe_get(qualifiers, "product")
                    name_parts = [self.abbreviation, self.sample_name, contig, locus_tag, gene, product,
                                  seq_feature.location]
                    out_name = ("|".join(Utils.remove_empty_values(
                        [re.sub(" +", "_", str(i).strip()) for i in name_parts])))
                    self._out_pfasta_records.append(
                        SeqRecord(Seq(Utils.safe_get(qualifiers, "translation"), IUPAC.protein), id=id_str,
                                  description=""))
                    self._out_annotations[id_str] = out_name

    def export(self):
        print("Save protein FASTA to: '{}'".format(self.out_pfasta))
        Utils.dump_string(s="".join([i.format("fasta") for i in self._out_pfasta_records]), file=self.out_pfasta)
        annotation_file = "{}_annotation.tsv".format(self.out_pfasta)
        print("Save protein annotation to: '{}'".format(annotation_file))
        Utils.dump_string(s="".join(["{}\t{}\n".format(
            k, self._out_annotations.get(k)) for k in self._out_annotations]), file=annotation_file)

    def convert(self):
        self.parse()
        self.export()


class Utils:
    @staticmethod
    def safe_findall(pattern, string, idx: int = 0):
        import re
        try:
            return re.findall(pattern, string)[idx]
        except IndexError:
            print("Warning! Can't find the regex pattern '{}' within the string: '{}'".format(pattern, string))
            return "unknown"

    @staticmethod
    def safe_get(d: dict, k: str):
        v = d.get(k)
        if not v:
            return "unknown"
        if isinstance(v, list):
            return "".join(v)
        return v

    @staticmethod
    def safe_extract_int(s: str):
        out = Utils.safe_findall("\d+", s)
        if len(out) > 0:
            return int(out)
        return 0

    @staticmethod
    def remove_empty_values(input_list):
        return [j for j in [i.strip() for i in input_list] if len(j) > 0]

    @staticmethod
    def dump_string(s: str, file: str):
        with open(file, mode="w", encoding="utf-8") as f:
            f.write(s)
            f.close()


if __name__ == '__main__':
    parser = ArgParser()
    converter = Converter()
    converter.convert()
