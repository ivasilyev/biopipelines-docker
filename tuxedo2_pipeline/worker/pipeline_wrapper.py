#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import argparse
import subprocess
import logging


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This script performs bowtie/bowtie2-based alignment and processes the output using MetaPhlAn 2. REFDATA files are generated by the 'cook_the_reference' script")
    starting_parser.add_argument("-s", "--sampledata", required=True,
                                 help="Input list containing two tab-delimited columns for colorspace or non-colorspace sequences and three for paired-end sequences: sample name and absolute path(s). May contain a header")
    starting_parser.add_argument("-h", "--hisat2_idx", required=True,
                                 help="HISAT2 reference index file mask")
    starting_parser.add_argument("-g", "--gtf", required=True,
                                 help="Annotated reference GTF file")
    starting_parser.add_argument("-m", "--mask", default="",
                                 help="(Optional) Mask to be added to resulting files")
    starting_parser.add_argument("-t", "--threads", default=None, type=int,
                                 help="(Optional) Number of CPU cores to use, maximal by default")
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


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        logging.info("Completed subprocess command: '" + str(input_direction_list) + "' with logging to STDOUT")
        return output.decode("utf-8").replace('\r', '').replace('\n', '')
    else:
        logging.info("Completed subprocess command: '" + str(input_direction_list) + "' with logging to '" + output_direction + "'")
        file_append(output.decode("utf-8"), output_direction)


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def file_to_list(file):
    file_buffer = open(file, 'rU')
    output_list = [j for j in [re.sub('[\r\n]', '', i) for i in file_buffer] if len(j) > 0]
    file_buffer.close()
    return output_list


def find_latest_changed_file(mask):
    return subprocess.getoutput("ls -1t -d " + mask + " | head -1")


def verify_path(path, file_or_dir):
    if path:
        if file_or_dir == "file" and os.path.isfile(path) is True:
            return os.path.abspath(path)
        if file_or_dir == "dir" and os.path.isdir(path) is True:
            return os.path.abspath(ends_with_slash(path))
        else:
            print("Cannot comprehend " + path + " as \"" + file_or_dir + "\"! Exiting...")
            sys.exit(2)


def parse_namespace():
    namespace = parse_args()
    for i in (namespace.sampledata, namespace.gtf):
        if not os.path.isfile(i):
            raise ValueError("Not found: '" + i + "'\nIf you're using Docker, please make sure you have mounted required volume with the '-v' flag.")
    default_threads = int(subprocess.getoutput("nproc"))
    if not namespace.threads or default_threads < namespace.threads:
        namespace.threads = default_threads
    is_path_exists(namespace.output)
    namespace.output = ends_with_slash(namespace.output)
    return namespace.sampledata, namespace.hisat2_idx, namespace.gtf, namespace.mask, str(namespace.threads), namespace.output


def create_dirs(dirs_list):
    for i in dirs_list:
        is_path_exists(i)


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


def per_sample_queue(sample_name, samples_paths_list):
    samples_paths_list = [i.strip() for i in samples_paths_list]
    for testing_path in samples_paths_list:
        if not os.path.isfile(testing_path):
            logging.warning("Cannot load the file: " + testing_path)
            samples_paths_list.remove(testing_path)
    if len(samples_paths_list) > 2:
        logging.warning("More than 2 paths supplied for sample '" + sample_name + "', using first two values of the list: ", str(samples_paths_list))
        samples_paths_list = samples_paths_list[:2]
    if len(samples_paths_list) < 2:
        logging.warning("Less than 2 paths supplied for sample '" + sample_name + "', cancelling processing the list: ", str(samples_paths_list))
        return
    # Now we have the verified list containing 2 NGS reads files
    output_dict = {sample_name: {}}
    output_dict[sample_name]["hisat2"] = outputDir + "hisat2/" + sample_name + "_" + inputMask + ".sam"
    external_route(["hisat2", "-p", cpuThreadsString, "--dta", "-x", hisat2IndexFileMask, "-1", samples_paths_list[0], "-2", samples_paths_list[1], "-S", output_dict[sample_name]["hisat2"]], outputDir + "logs/" + sample_name + "_" + inputMask + "_hisat2.log")
    output_dict[sample_name]["samtools"] = outputDir + "samtools/" + sample_name + "_" + inputMask + ".bam"
    external_route(["samtools", "sort", "-@", cpuThreadsString, "-o", output_dict[sample_name]["samtools"],  output_dict[sample_name]["hisat2"]], outputDir + "logs/" + sample_name + "_" + inputMask + "_samtools.log")
    output_dict[sample_name]["stringtie"] = outputDir + "stringtie/" + sample_name + "_" + inputMask + ".gtf"
    external_route(["stringtie", "-p", cpuThreadsString, referenceGTFFileName, "-o", output_dict[sample_name]["stringtie"], "-l", output_dict[sample_name]["samtools"]], outputDir + "logs/" + sample_name + "_" + inputMask + "_stringtie.log")
    file_append(output_dict[sample_name]["stringtie"], outputDir + "stringtie/" + filename_only(sampleDataFileName) + "_" + inputMask + "_gtfs_list.txt")
    return output_dict


def merging_queue(input_dict):
    merged_gtfs_file_name = outputDir + "stringtie/" + filename_only(sampleDataFileName) + "_" + inputMask + "_merged.gtf"
    external_route(["stringtie", "--merge", "-p", cpuThreadsString, "-G", referenceGTFFileName, "-o", merged_gtfs_file_name, outputDir + "stringtie/" + filename_only(sampleDataFileName) + "_" + inputMask + "_gtfs_list.txt"], outputDir + "logs/" + filename_only(sampleDataFileName) + "_" + inputMask + "_stringtie-merge.log")
    external_route(["gffcompare", "-r", referenceGTFFileName, "-G", "-o", "prefix", merged_gtfs_file_name], outputDir + "logs/" + filename_only(sampleDataFileName) + "_" + inputMask + "_gffcompare.log")
    for sample_name in input_dict:
        external_route(["stringtie", "-e", "-B", "-p", cpuThreadsString, "-G", merged_gtfs_file_name, "-o", outputDir + "ballgown/" + sample_name + "_" + inputMask + ".gtf", input_dict[sample_name]["samtools"]], outputDir + "logs/" + sample_name + "_" + inputMask + "_stringtie-e.log")


if __name__ == '__main__':
    sampleDataFileName, hisat2IndexFileMask, referenceGTFFileName, inputMask, cpuThreadsString, outputDir = parse_namespace()
    scriptDir = ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))
    create_dirs([outputDir + i for i in ["hisat2", "samtools", "stringtie", "gffcompare", "ballgown", "logs"]])
    currentTime = get_time()
    logging.basicConfig(format=u'%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG, filename=outputDir + "logs/" + "_".join([filename_only(sys.argv[0]), inputMask, currentTime, "master.log"]))
    logging.info("Command: " + ' '.join(str(i) for i in sys.argv))
    logging.info("Main process ID: " + str(os.getpid()))
    logging.info("Using CPU threads: " + cpuThreadsString)
    queueDict = {}
    for singleSampleDataLine in file_to_list(sampleDataFileName):
        try:
            sampleName, samplePathsList = singleSampleDataLine.split("\t")[0], singleSampleDataLine.split("\t")[1:]
            singleSampleDict = per_sample_queue(sampleName, samplePathsList)
            if isinstance(singleSampleDict, dict):
                queueDict.update(singleSampleDict)
        except IndexError:
            logging.warning("Cannot parse the sampledata string: " + singleSampleDataLine)
            continue
    logging.info("Number of successfully processed samples: " + str(len(queueDict)))
    merging_queue(queueDict)
    logging.info("COMPLETED")
    print("COMPLETED")