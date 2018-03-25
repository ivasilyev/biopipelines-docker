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


def parse_namespace():
    namespace = parse_args()
    if os.path.isfile(namespace.cfg):
        with open(namespace.cfg, 'r') as stream:
            parsed_cfg = yaml.load(stream)
    else:
        if namespace.cfg.lower().startswith("http"):
            parsed_cfg = yaml.load(external_route("curl", "-fsSL", namespace.cfg))
        else:
            raise ValueError("Cannot parse config file or URL!")
    is_path_exists(namespace.output)
    return parsed_cfg, ends_with_slash(namespace.output)


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    return output.decode("utf-8")


def template2yaml(template_url, output_file_name):
    template = jinja2.Template(template_url)
    with open(output_file_name, 'w') as stream:
        yaml.dump(template.render(yaml.load(parsedCFG)), stream, default_flow_style=False)


if __name__ == '__main__':
    parsedCFG, outputDir = parse_namespace()
    # Paste template URLs here
    exportDict = {"master": "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage-master.yaml",
                  "worker": "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage-worker.yaml"}
    for templateName in exportDict:
        template2yaml(exportDict[templateName], outputDir + templateName + ".yaml")
