"""
Find kmers of certain length that repeat in genome
Determine if they appear in Illumina assembly or break regions
"""
sample = config['genomes']
mersize = config['mersize']
mercount = config['mercount']
nmag_bed = config['nmag_bed']
imag_bed = config['imag_bed']
prefix = sample.replace(".fa", "")

rule all:
    input:
        expand("{prefix}.mer_counts{mersize}.jf", prefix = prefix, mersize = mersize),
        expand("{prefix}.mer_counts{mersize}_dumps.fa", prefix = prefix, mersize = mersize),
        expand("{prefix}.mer_counts{mersize}.histo", prefix = prefix, mersize = mersize),
        expand("{prefix}.{mercount}repeats{mersize}.fa", prefix = prefix, mercount = mercount, mersize = mersize),
        expand("{prefix}.fa.sa", prefix = prefix),
        expand("{prefix}.{mercount}repeats{mersize}_alignment.bam.bai", prefix = prefix, mercount = mercount, mersize = mersize),
        #expand("{prefix}_uncovered.fasta", prefix = prefix),
        expand("{prefix}_unaligned.fasta", prefix = prefix)

# find counts of all kmers of a given size
rule kmer_counts:
    input:
        #fa = {sample}
        "{prefix}.fa"
    output:
        "{prefix}.mer_counts{mersize}.jf"
    shell:
        # "jellyfish count -m {mersize} -s 1M -t 10 -C -o \
        # {prefix}.mer_counts{mersize}.jf {sample}"
        "jellyfish count -m {mersize} -s 1M -t 10 -C -o \
        {prefix}.mer_counts{mersize}.jf {prefix}.fa"

# output all kmers to fasta file with associated counts
rule counts_to_fasta:
    input:
        "{prefix}.mer_counts{mersize}.jf"
    output:
        "{prefix}.mer_counts{mersize}_dumps.fa"
    shell:
        "jellyfish dump {prefix}.mer_counts{mersize}.jf > \
        {prefix}.mer_counts{mersize}_dumps.fa"

# generate histogram of kmer counts
rule hist_mer_counts:
    input:
        "{prefix}.mer_counts{mersize}.jf"
    output:
        "{prefix}.mer_counts{mersize}.histo"
    shell:
        "jellyfish histo {prefix}.mer_counts{mersize}.jf > \
        {prefix}.mer_counts{mersize}.histo"

# pull all kmers of a given size that occur with at least a given frequency
rule pull_mers:
    input:
        "{prefix}.mer_counts{mersize}_dumps.fa"
    output:
        "{prefix}.{mercount}repeats{mersize}.fa"
    shell:
        "python scripts/pull_mers.py {prefix}.mer_counts{mersize}_dumps.fa \
        {prefix}.{mercount}repeats{mersize}.fa {mercount}"

rule index_genome:
    input:
        fa = {sample}
    output:
        "{prefix}.fa.amb",
        "{prefix}.fa.ann",
        "{prefix}.fa.bwt",
        "{prefix}.fa.pac",
        "{prefix}.fa.sa"
    shell:
        #"bwa index {sample}"
        "bwa index {prefix}.fa"

# align reads back to genome sequence
rule align_kmers:
    input:
        "{prefix}.fa.amb",
        "{prefix}.fa.ann",
        "{prefix}.fa.bwt",
        "{prefix}.fa.pac",
        "{prefix}.fa.sa",
        "{prefix}.{mercount}repeats{mersize}.fa"
    output:
        "{prefix}.{mercount}repeats{mersize}_alignment.bam",
        "{prefix}.{mercount}repeats{mersize}_alignment.bam.bai"
    threads: 12
    resources:
        time = 1,
        mem = 24
    shell:
        # "bwa mem -a -t 12 {sample} {prefix}.{mercount}repeats{mersize}.fa | \
        # samtools sort --thread 24 > {prefix}.{mercount}repeats{mersize}_alignment.bam; \
        # samtools index {prefix}.{mercount}repeats{mersize}_alignment.bam \
        # {prefix}.{mercount}repeats{mersize}_alignment.bam.bai"
        "bwa mem -a -t 12 {prefix}.fa {prefix}.{mercount}repeats{mersize}.fa | \
        samtools sort --thread 24 > {prefix}.{mercount}repeats{mersize}_alignment.bam; \
        samtools index {prefix}.{mercount}repeats{mersize}_alignment.bam \
        {prefix}.{mercount}repeats{mersize}_alignment.bam.bai"

# get bed file of regions that are not covered by illumina assembly
rule unaligned_nmag_bed:
    input:
        nmag_bed = {config['nmag_bed']},
        imag_bed = {config['imag_bed']}
    output:
        "{prefix}_uncovered.bed"
    shell:
        "bedtools subtract -a {input.nmag_bed} -b {input.imag_bed} > {prefix}_uncovered.bed"

# FIXME: breaks at this rule - RuleException TypeError expected str got bool
# get corresponding fasta of regions that are not covered by illumina assembly
rule unaligned_nmag_fasta:
    input:
        "{prefix}_uncovered.bed"
    output:
        "{prefix}_unaligned.fasta"
    shell:
        #"bedtools getfasta -fi {sample} -bed {prefix}_uncovered.bed | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > {prefix}_unaligned.fasta"
        "bedtools getfasta -fi {prefix}.fa -bed {prefix}_uncovered.bed | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > {prefix}_unaligned.fasta"
