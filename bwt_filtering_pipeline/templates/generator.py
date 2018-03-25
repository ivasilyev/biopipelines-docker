#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import yaml
import jinja2
import os


def parse_args():
    starting_parser = argparse.ArgumentParser(description="The generator for 'bwt_filtering_pipeline' templates")
    starting_parser.add_argument("-c", "--cfg", required=True,
                                 help="Configuration file or URL")
    starting_parser.add_argument("-m", "--master", required=True,
                                 help="MASTER file or URL")
    starting_parser.add_argument("-w", "--worker", required=True,
                                 help="WORKER file or URL")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory")
    return starting_parser.parse_args()


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def read_file_or_url(path):
    if os.path.isfile(path):
        with open(path, 'r') as file:
            return file.read()
    else:
        if path.startswith("http") or path.startswith("www."):
            return external_route("curl", "-fsSL", path)
        else:
            raise ValueError("Cannot load  file or URL!")


def parse_namespace():
    namespace = parse_args()
    is_path_exists(namespace.output)
    namespace.cfg, namespace.master, namespace.worker = [read_file_or_url(i) for i in [namespace.cfg, namespace.master, namespace.worker]]
    return namespace.cfg, namespace.master, namespace.worker, ends_with_slash(namespace.output)


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    return output.decode("utf-8")


def template2yaml(buffered_template, output_file_name):
    template = jinja2.Template(buffered_template)
    with open(output_file_name, 'w') as stream:
        yaml.dump(template.render(yaml.load(bufferedCFG)), stream, default_flow_style=False)


if __name__ == '__main__':
    bufferedCFG, bufferedMaster, bufferedWorker, outputDir = parse_namespace()
    # Paste template URLs here
    exportDict = {"master": bufferedMaster,
                  "worker": bufferedWorker}
    for templateName in exportDict:
        template2yaml(exportDict[templateName], outputDir + templateName + ".yaml")
