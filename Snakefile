#!/usr/bin/env python3

import os
import pathlib
import pandas


#############
# FUNCTIONS #
#############

def resolve_path(my_path):
    return str(pathlib.Path(my_path).resolve())


# def find_input_files(wildcards):
def find_input_files(wildcards):
    my_treatment = wildcards.treatment
    my_rep = wildcards.rep
    my_samples = list(
        sample_key[(sample_key.treatment == my_treatment) &
                   (sample_key.replicate == int(my_rep))]['file_name'])
    my_r1 = [os.path.join(read_dir, x)
             for x in my_samples if '_R1_' in x][0]
    my_r2 = [os.path.join(read_dir, x)
             for x in my_samples if '_R2_' in x][0]
    return {'r1': my_r1, 'r2': my_r2}

###########
# GLOBALS #
###########

sample_key_file = 'data/sample_key.csv'
read_dir = 'data/fastq'
bbduk_adaptors = 'venv/bin/resources/adapters.fa'
bbduk_contaminants = 'venv/bin/resources/sequencing_artifacts.fa.gz'

#########
# SETUP #
#########

# generate name to filename dictionary
sample_key = pandas.read_csv(sample_key_file)

#########
# RULES #
#########

rule target:
    input:
        expand('output/trim_clip/{treatment}_{rep}_r{r}.fastq',
               treatment=['trt1', 'trt2', 'untreated'],
               rep=['1', '2'],
               r=['1', '2'])

# 1. trim and clip with bbduk
rule trim_clip:
    input:
        unpack(find_input_files),
        adaptors = bbduk_adaptors,
        contaminants = bbduk_contaminants
    output:
        r1 = 'output/trim_clip/{treatment}_{rep}_r1.fastq',
        r2 = 'output/trim_clip/{treatment}_{rep}_r2.fastq'
    log:
        trim_log = 'output/logs/{treatment}_{rep}_trim.log',
        trim_stats = 'output/trim_clip/{treatment}_{rep}_trim-stats.txt',
        filter_log = 'output/logs/{treatment}_{rep}_filter.log',
        filter_stats = 'output/trim_clip/{treatment}_{rep}_filter-stats.txt'
    threads:
        10
    shell:
        'bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in=\'{input.r1}\' '
        'in2=\'{input.r2}\' '
        'out=stdout.fastq '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'ref={input.adaptors} '
        'stats={log.trim_stats} '
        'statscolumns=5 '
        '2> {log.trim_log} '
        '| bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in=stdin.fastq '
        'interleaved=t '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={input.contaminants} '
        'k=31 hdist=1 stats={log.filter_stats} '
        'minlength=50 '
        '2> {log.filter_log}'

