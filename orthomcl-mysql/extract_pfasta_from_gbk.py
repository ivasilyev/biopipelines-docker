#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
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
        _parser.add_argument("-s", "--sample_name", metavar="<str>", type=str, default="",
                             help="Sample name to add into prefix")
        _parser.add_argument("-o", "--output", metavar="<output.faa>", type=str, required=True, help="Output file name")
        self._namespace = _parser.parse_args()
        self.input_gbk = self._namespace.input
        self.sample_name = self._namespace.sample_name
        self.out_pfasta = self._namespace.output


class Converter:
    def __init__(self):
        self.input_gbk = parser.input_gbk
        self.sample_name = parser.sample_name
        self.out_pfasta = parser.out_pfasta
        self._out_pfasta_records = []
        self._out_annotations = pd.DataFrame()

    def parse(self):
        seq_records = list(SeqIO.parse(self.input_gbk, "genbank"))
        max_ids = 0
        for seq_record in seq_records:
            for seq_feature in seq_record.features:
                if seq_feature.type == "CDS" and "translation" in seq_feature.qualifiers.keys():
                    max_ids += 1
        _id = 0
        annotations = []
        # LOCUS
        for seq_record in seq_records:
            # CDS
            for seq_feature in seq_record.features:
                qualifiers = seq_feature.qualifiers.copy()
                if seq_feature.type == "CDS" and "translation" in qualifiers.keys():
                    _id += 1
                    id_str = "PFASTA_ID_{}".format(str(_id).zfill(len(str(max_ids))))
                    contig = "contig_{}".format(seq_record.id)
                    locus_tag = "CDS_{}".format(Utils.safe_get(qualifiers, "locus_tag"))
                    gene = Utils.safe_get(qualifiers, "gene")
                    product = Utils.safe_get(qualifiers, "product")
                    pfasta_seq = Seq(Utils.safe_get(qualifiers, "translation"), IUPAC.protein)
                    annotation = {"sample_name": self.sample_name, "contig": contig, "locus_tag": locus_tag,
                                  "gene": gene, "product": product, "location": str(seq_feature.location),
                                  "pfasta_id": id_str, "nfasta_bp": len(seq_feature), "pfasta_bp": len(pfasta_seq)}
                    self._out_pfasta_records.append(SeqRecord(pfasta_seq, id=id_str, description=""))
                    annotations.append(annotation)
        self._out_annotations = pd.DataFrame(annotations)

    def export(self):
        print("Save protein FASTA to: '{}'".format(self.out_pfasta))
        Utils.dump_string(s="".join([i.format("fasta") for i in self._out_pfasta_records]), file=self.out_pfasta)
        annotation_file = "{}_annotation.tsv".format(".".join(self.out_pfasta.split(".")[:-1]))
        print("Save protein annotation to: '{}'".format(annotation_file))
        self._out_annotations.to_csv(annotation_file, encoding="utf-8", sep="\t", index=False, header=True)

    def convert(self):
        self.parse()
        self.export()


class Utils:
    @staticmethod
    def safe_findall(pattern, string, idx: int = 0):
        try:
            return re.findall(pattern, string)[idx]
        except IndexError:
            print("Warning! Can't find the regex pattern '{}' within the string: '{}'".format(pattern, string))
            return ""

    @staticmethod
    def safe_get(d: dict, k: str):
        v = d.get(k)
        if not v:
            return ""
        if isinstance(v, list):
            return "".join(v)
        return v

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
