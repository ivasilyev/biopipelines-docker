#!/usr/bin/env python
# The original topic: https://bitbucket.org/jason_c_kwan/orthomcl_docker/src/master/

import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Run OrthoMCL on a group of proteomes.',
                                 epilog='Stages: 1 - Preprocessing (orthomclAdjustFasta and orthomclFilterFasta), 2 - BLAST search, 3 - Database processing, 4 - Post processing (MCL)')
parser.add_argument('-t', '--table_path', metavar='<input_table.tab>',
                    help='Path to a table containing the relative paths of proteome fastas (column 1) and species abbreviations to use (column 2)',
                    required=True)
parser.add_argument('-s', '--start_stage', metavar='<1|2|3|4>', help='Stage to start the pipeline', type=int,
                    choices=[1, 2, 3, 4], default=1)
parser.add_argument('-e', '--end_stage', metavar='<1|2|3|4>', help='Stage to end the pipeline', type=int,
                    choices=[1, 2, 3, 4], default=4)
parser.add_argument('-p', '--processors', metavar='<int>',
                    help='Number of processors to use (only really relevant for BLAST stage)', type=int, default=1)
parser.add_argument('-c', '--config_file', metavar='<file>', help='Path to config file to use with orthomcl',
                    required=True)

args = vars(parser.parse_args())
table_path = os.path.abspath(args['table_path'])
# Note, assumes the output dir is the one that holds the table
output_dir = '/'.join(table_path.split('/')[:-1])
os.chdir(output_dir)
start_stage = args['start_stage']
end_stage = args['end_stage']
num_threads = args['processors']
config_file = os.path.abspath(args['config_file'])

# Some error checks
if start_stage > end_stage:
    print('Error! Start stage must be before end stage')
    exit(1)

stages_to_do = range(start_stage, end_stage + 1)

if 1 in stages_to_do:
    # Parse input table (format proteome_path\tabbreviation)
    abbreviations = dict()
    with open(table_path) as table:
        for line in table:
            line_list = line.rstrip().split('\t')
            proteome_path = os.path.abspath(line_list[0])
            abbreviation = line_list[1]
            abbreviations[proteome_path] = abbreviation

    # orthomclAdjustFasta
    os.system('mkdir compliantFasta')
    os.chdir('compliantFasta')
    for proteome in abbreviations:
        command = 'orthomclAdjustFasta ' + abbreviations[proteome] + ' ' + proteome + ' 1'
        print(command)
        os.system(command)
    os.chdir('../')

    # orthomclFilterFasta
    os.system('orthomclFilterFasta compliantFasta 10 20')

if 2 in stages_to_do:
    # Make diamond database
    command = 'diamond makedb --in goodProteins.fasta -d goodProteins.fasta'
    print(command)
    os.system(command)

    # Do diamond search
    command = 'diamond blastp --threads ' + str(
        num_threads) + ' --db goodProteins.fasta --outfmt 6 --out all_v_all.blastp --query goodProteins.fasta --max-target-seqs 100000 --evalue 1e-5 --masking 1'
    print(command)
    os.system(command)

    # Process blast results
    command = 'orthomclBlastParser all_v_all.blastp compliantFasta >> similarSequences.txt'
    print(command)
    os.system(command)

if 3 in stages_to_do:
    # Set up database
    command = 'mysql -u root -pPAssw0rd -e "create database orthomcl"'
    print(command)
    os.system(command)
    command = 'orthomclInstallSchema ' + config_file + ' orthomclInstallSchema_sql.log'
    print(command)
    os.system(command)

    # Load data into database
    command = 'orthomclLoadBlast ' + config_file + ' similarSequences.txt'
    print(command)
    os.system(command)

    # Process pairs
    command = 'orthomclPairs ' + config_file + ' orthomclPairs.log cleanup=yes'
    print(command)
    os.system(command)

    # Dump output
    command = 'orthomclDumpPairsFiles ' + config_file
    print(command)
    os.system(command)

if 4 in stages_to_do:
    # Run MCL
    command = 'mcl mclInput --abc -I 1.5 -o mclOutput'
    print(command)
    os.system(command)

    # Make groups file
    command = 'orthomclMclToGroups groups 1000 < mclOutput > groups.txt'
    print(command)
    os.system(command)
