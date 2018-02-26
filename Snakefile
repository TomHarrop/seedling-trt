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
star_reference_folder = 'output/star_reference'

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
        expand(('output/star_pass2/'
                '{treatment}_{rep}.Aligned.sortedByCoord.out.bam.bai'),
               treatment=['trt1', 'trt2', 'untreated'],
               rep=['1', '2'])

# 4. plots
rule plot_counts_stats:
    input:
        bam_files = expand(('output/star_pass2/'
                '{treatment}_{rep}.Aligned.sortedByCoord.out.bam.bai'),
               treatment=['trt1', 'trt2', 'untreated'],
               rep=['1', '2']),
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        counts_plot = 'output/mapping_stats/counts_per_category.pdf',
        intron_exon_plot = 'output/mapping_stats/intron_exon_counts.pdf',
        feature_counts = 'output/mapping_stats/feature_counts.csv'
    threads:
        6
    script:
        'src/count_reads_per_feature.R'


# 3. map
rule index_bam_files:
    input:
        ('output/star_pass2/'
         '{treatment}_{rep}.Aligned.sortedByCoord.out.bam')
    output:
        ('output/star_pass2/'
         '{treatment}_{rep}.Aligned.sortedByCoord.out.bam.bai')
    threads:
        1
    shell:
        'samtools index {input} {output}'


rule star_second_pass:
    input:
        r1 = 'output/trim_clip/{treatment}_{rep}_r1.fastq',
        r2 = 'output/trim_clip/{treatment}_{rep}_r2.fastq',
        star_reference = 'output/star_reference/Genome',
        junctions = expand(
            expand('output/star_pass1/{treatment}_{rep}.SJ.out.tab',
                   treatment=['trt1', 'trt2', 'untreated'],
                   rep=['1', '2']))
    output:
        bam = ('output/star_pass2/'
               '{treatment}_{rep}.Aligned.sortedByCoord.out.bam'),
        counts = 'output/star_pass2/{treatment}_{rep}.ReadsPerGene.out.tab'
    threads:
        15
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star_pass2/{treatment}_{rep}.'
    log:
        'output/logs/STAR_pass2_{treatment}_{rep}.log'
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
        r1 = 'output/trim_clip/{treatment}_{rep}_r1.fastq',
        r2 = 'output/trim_clip/{treatment}_{rep}_r2.fastq',
        star_reference = 'output/star_reference/Genome'
    output:
        sjdb = 'output/star_pass1/{treatment}_{rep}.SJ.out.tab'
    threads:
        15
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star_pass1/{treatment}_{rep}.'
    log:
        'output/logs/STAR_pass1_{treatment}_{rep}.log'
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
        'in={input.r1} '
        'in2={input.r2} '
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
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 149 '
        '&> {log}'


