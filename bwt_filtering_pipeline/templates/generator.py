#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import yaml
import jinja2
import os


class Initializer:
    def __init__(self):
        self._namespace = self.parse_args()
        self.config = self._namespace.config
        self.master = self._namespace.master
        self.worker = self._namespace.worker
        self.output = self._namespace.output
    @staticmethod
    def parse_args():
        starting_parser = argparse.ArgumentParser(description="The generator for 'bwt_filtering_pipeline' templates")
        starting_parser.add_argument("-c", "--config", required=True,
                                     help="Configuration file or URL")
        starting_parser.add_argument("-m", "--master", required=True,
                                     help="MASTER file or URL")
        starting_parser.add_argument("-w", "--worker", required=True,
                                     help="WORKER file or URL")
        starting_parser.add_argument("-o", "--output", required=True,
                                     help="Output directory")
        return starting_parser.parse_args()


class KubernetesJobChartsGenerator:
    def __init__(self, config, master, worker):
        self._bufferedCFG = self._read_file_or_url(config)
        self._exportDict = {"master": self._read_file_or_url(master),
                            "worker": self._read_file_or_url(worker)}
    @staticmethod
    def _read_file_or_url(path):
        if os.path.isfile(path):
            open(file=path, mode="r", encoding="utf-8")
            with open(path, 'r') as f:
                return f.read()
        else:
            if path.startswith("http") or path.startswith("www."):
                return subprocess.getoutput("curl -fsSL {}".format(path))
            else:
                raise ValueError("Cannot load  file or URL: '{}'".format(path))
    def _template2yaml(self, buffered_template):
        template = jinja2.Template(buffered_template)
        return "---\n{}\n".format(template.render(yaml.load(self._bufferedCFG)))
    def export_yaml(self, output_dir):
        output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]
        for template_name in self._exportDict:
            s = self._template2yaml(self._exportDict[template_name])
            with open(file="{a}{b}.yaml".format(a=output_dir, b=template_name), mode="w", encoding="utf-8") as f:
                f.write(s)


if __name__ == '__main__':
    mainInitializer = Initializer()
    generator = KubernetesJobChartsGenerator(config=mainInitializer.config,
                                             master=mainInitializer.master,
                                             worker=mainInitializer.worker)
    generator.export_yaml(mainInitializer.output)
