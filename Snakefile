#!/usr/bin/env python3

import os
import pathlib
import pandas
import shlex


#############
# FUNCTIONS #
#############

def resolve_path(my_path):
    return str(pathlib.Path(my_path).resolve())

# find acceptable combos from sample_key
def find_combos(sample_key, experiment):
    my_key = sample_key[sample_key['experiment'] == experiment]
    my_samples = sorted(set(my_key['sample']))
    my_replicates = sorted(set(my_key['replicate']))
    return({'samples': my_samples,
            'replicates': my_replicates})

# def find_input_files(wildcards):
def find_input_files(wildcards):
    my_experiment = wildcards.experiment
    my_sample = wildcards.sample
    my_replicate = wildcards.replicate
    my_key = sample_key[
        (sample_key['experiment'] == my_experiment) &
        (sample_key['sample'] == my_sample) &
        (sample_key['replicate'] == my_replicate)]
    my_r1 = my_key.iloc[0]['r1_rel_path']
    my_r2 = my_key.iloc[0]['r2_rel_path']
    return {'r1': my_r1, 'r2': my_r2}


###########
# GLOBALS #
###########

sample_key_file = 'data/sample_key.csv'
read_dir = 'data/fastq'
bbduk_adaptors = '/adapters.fa'                     # in the bbmap container
bbduk_contaminants = '/sequencing_artifacts.fa.gz'  # in the bbmap container
star_reference_folder = 'output/star_reference'

# containers
bbduk_container = ('shub://TomHarrop/singularity-containers:bbmap_38.00'
                   '@a773baa8cc025cc5b5cbee20e507fef7')
star_container = ('shub://TomHarrop/singularity-containers:star_2.6.0c'
                  '@eaa90a258fdb26b6b0ce7d07246ffe2c')
bioc_container = ('shub://TomHarrop/singularity-containers:bioconductor_3.7'
                  '@2785c89cc4ef1cbb08c06dde9ecf9544')



#########
# SETUP #
#########

# read the sample key
sample_key = pandas.read_csv(sample_key_file)

# get all experiments
all_expts = sorted(set(sample_key['experiment']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/mapping_stats/{experiment}/feature_counts.csv',
               experiment=all_expts)

# # 4. plots
rule plot_counts_stats:
    input:
        bam_files = lambda wildcards: expand(
            ('output/star_pass2/{0}/{{sample}}_{{replicate}}.'
             'Aligned.sortedByCoord.out.bam').format(wildcards.experiment),
            sample=find_combos(sample_key, wildcards.experiment)['samples'],
            replicate=find_combos(sample_key, wildcards.experiment)['replicates']),
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        feature_counts = 'output/mapping_stats/{experiment}/feature_counts.csv'
    log:
        log = 'output/logs/plot_counts_stats_{experiment}.log'
    threads:
        10
    singularity:
        bioc_container
    script:
        'src/count_reads_per_feature.R'


# # 3. map
rule star_second_pass:
    input:
        r1 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r1.fastq',
        r2 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r2.fastq',
        star_reference = 'output/star_reference/Genome',
        junctions = lambda wildcards: expand(
            ('output/star_pass1/{0}/{{sample}}_{{replicate}}.'
             'SJ.out.tab').format(wildcards.experiment),
            sample=find_combos(sample_key, wildcards.experiment)['samples'],
            replicate=find_combos(sample_key, wildcards.experiment)['replicates'])
    output:
        bam = 'output/star_pass2/{experiment}/{sample}_{replicate}.Aligned.sortedByCoord.out.bam',
        counts = 'output/star_pass2/{experiment}/{sample}_{replicate}.ReadsPerGene.out.tab'
    threads:
        30
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star_pass2/{experiment}/{sample}_{replicate}.'
    log:
        'output/logs/STAR_pass2_{experiment}_{sample}_{replicate}.log'
    priority:
        1
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--outReadsUnmapped Fastx '
        '--quantMode GeneCounts '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_first_pass:
    input:
        r1 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r1.fastq',
        r2 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r2.fastq',
        star_reference = 'output/star_reference/Genome'
    output:
        sjdb = 'output/star_pass1/{experiment}/{sample}_{replicate}.SJ.out.tab'
    threads:
        30
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star_pass1/{experiment}/{sample}_{replicate}.'
    log:
        'output/logs/STAR_pass1/{experiment}_{sample}_{replicate}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

# 2. trim and clip with bbduk
rule trim_clip:
    input:
        unpack(find_input_files),
    output:
        r1 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r1.fastq',
        r2 = 'output/trim_clip/{experiment}/{sample}_{replicate}_r2.fastq'
    params:
        adaptors = bbduk_adaptors,
        contaminants = bbduk_contaminants
    log:
        trim_log = 'output/logs/{experiment}_{sample}_{replicate}_trim.log',
        trim_stats = 'output/trim_clip/{experiment}/{sample}_{replicate}_trim-stats.txt',
        filter_log = 'output/logs/{experiment}_{sample}_{replicate}_filter.log',
        filter_stats = 'output/trim_clip/{experiment}/{sample}_{replicate}_filter-stats.txt'
    threads:
        2
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'ref={params.adaptors} '
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
        'ref={params.contaminants} '
        'k=31 hdist=1 stats={log.filter_stats} '
        'minlength=50 '
        '2> {log.filter_log}'

# 1. process reference for STAR
rule star_reference:
    input:
        fasta = 'data/ref/TAIR10_Chr.all.fasta',
        gtf = 'data/ref/Araport11_GFF3_genes_transposons.201606.gtf'
    output:
        'output/star_reference/Genome'
    params:
        genome_dir = star_reference_folder
    threads:
        30
    log:
        'output/logs/star_reference.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 149 '
        '&> {log}'


