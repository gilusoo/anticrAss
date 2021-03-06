'''

Cross assembly pipeline for bacteriophage discovery from metagenomic samples.

To run this pipeline on antihll, use:

	snakemake -s anticrass.snakefile --cluster "qsub -cwd -V" --latency-wait 60

'''
from os.path import join
from os import listdir
import os.path


configfile: "anticrass.yaml"

IDS = []
with open("SraAccList.txt", 'r') as infile:
	for line in infile:
		IDS.append(line.strip())

FASTQDIR = config["directories"]["fastq_files"]
SPADESDIR = config["directories"]["assembly"]
BOWTIEDIR = config["directories"]["bowtie"]
BLASTDIR = config["directories"]["blast_results"]
BLASTDB = config["databases"]["blast"]
CONCOCTDIR = config["directories"]["concoct_output"]
BENCHDIR = config["directories"]["benchmark"]

SPADES = config["commands"]["spades"]
BLASTREADER = config["commands"]["blast_reader"]

rule all:
	input: 
		expand(join(FASTQDIR, "{sample}_pass_1.fastq"), sample=IDS), 
		join(FASTQDIR, "AllReads_1.fastq"),
		join(FASTQDIR, "AllReads_2.fastq"),
		join(SPADESDIR, "contigs.fasta"),
		expand(join(BOWTIEDIR, "contigs_idx.{index}.bt2"), index=range(1,5)),
		expand(join(BOWTIEDIR, "contigs_idx.rev.{index}.bt2"), index=range(1,3)),
		expand(join(BOWTIEDIR, "SAM_files/{sample}.sam"), sample=IDS),
		expand(join(BOWTIEDIR, "BAM_files", "{sample}.bam"), sample=IDS),
		expand(join(BOWTIEDIR, "BAM_files/{sample}.bam.bai"), sample=IDS),
		join(CONCOCTDIR, "contigs_1k.bed"),
                join(CONCOCTDIR, "contigs_1k.fa"),
		join(CONCOCTDIR, "coverage_table.tsv"),
		join(CONCOCTDIR, "clustering_gt1000.csv"),
		join(CONCOCTDIR, "clustering_merged.csv"),
		expand(join(CONCOCTDIR, "fasta_bins/{bin}.fa"), bin=range(0,8)),
		expand(join(BLASTDIR, "{bin}.out"), bin=range(0,8)),
		expand(join(BLASTDIR, "protein_hits/{bin}_blast.tsv"), bin=range(0,8))

rule download:
	output:
		join(FASTQDIR, "{sample}_pass_1.fastq"),
		join(FASTQDIR, "{sample}_pass_2.fastq")
	params:
		sra_acc = "{sample}",
		outdir = FASTQDIR
	benchmark:
		join(BENCHDIR, "{sample}.download.benchmark.txt")
	shell:
		"fastq-dump --outdir {params.outdir} -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {params.sra_acc}"

rule concat1:
	input:
		expand(join(FASTQDIR, "{sample}_pass_1.fastq"), sample=IDS)
	output:
		join(FASTQDIR, "AllReads_1.fastq")
	benchmark:
                join(BENCHDIR, "concat1.benchmark.txt")
	shell:
		"cat {input} > {output}"

rule concat2:
        input:
                expand(join(FASTQDIR, "{sample}_pass_2.fastq"), sample=IDS)
        output:
                join(FASTQDIR, "AllReads_2.fastq")
	benchmark:
		join(BENCHDIR, "concat2.benchmark.txt")
	shell:
                "cat {input} > {output}"

rule assemble:
	input:
		join(FASTQDIR, "AllReads_1.fastq"),
		join(FASTQDIR, "AllReads_2.fastq")
	output:
		join(SPADESDIR, "contigs.fasta"),
		join(SPADESDIR, "dataset.info"), 
		join(SPADESDIR, "spades.log")
	params:
		spades_assemble = SPADES,
		outdir = config["directories"]["assembly"]
	benchmark:
		join(BENCHDIR, "assemble.benchmark.txt")
	shell:
		"python {params.spades_assemble} -m 1500 -1 {input[0]} -2 {input[1]} --threads 16 --only-assembler -o {params.outdir}"

rule bowtie_idx:
	input:
		join(SPADESDIR, "contigs.fasta")
	output:
		expand(join(BOWTIEDIR, "contigs_idx.{index}.bt2"), index=range(1,5)),
		expand(join(BOWTIEDIR, "contigs_idx.rev.{index}.bt2"), index=range(1,3))	
	params:
		idx = join(BOWTIEDIR, "contigs_idx")
	benchmark:
		join(BENCHDIR, "bowtieidx.benchmark.txt")
	shell:
		"bowtie2-build {input} {params.idx}" 

rule bowtie_scan:
	input:
		expand(join(BOWTIEDIR, "contigs_idx.{index}.bt2"), index=range(1,5)),
		expand(join(BOWTIEDIR, "contigs_idx.rev.{index}.bt2"), index=range(1,3))
	params:
		index = join(BOWTIEDIR, "contigs_idx"),
		forward = join(FASTQDIR, "{sample}_pass_1.fastq"),
		backward = join(FASTQDIR, "{sample}_pass_2.fastq")
	output:
		join(BOWTIEDIR, "SAM_files/{sample}.sam")
	benchmark:
		join(BENCHDIR, "{sample}.bowtiescan.benchmark.txt")
	shell:
		"bowtie2 -q -x {params.index} -1 {params.forward} -2 {params.backward} -S {output}"

rule sam_to_bam:
	input:
		join(BOWTIEDIR, "SAM_files/{sample}.sam")
	output:
		join(BOWTIEDIR, "BAM_files/{sample}.bam")
	benchmark:
		join(BENCHDIR, "{sample}.samtobam.benchmark.txt")
	shell:
		"samtools view -bS {input} | samtools sort -o {output}"

rule index_bam:
	input:
		join(BOWTIEDIR, "BAM_files/{sample}.bam")
	output:
		join(BOWTIEDIR, "BAM_files/{sample}.bam.bai")
	benchmark:
		join(BENCHDIR, "{sample}.indexbam.benchmark.txt")
	shell:
		"samtools index {input}"

'''
rule check_bam:
	input:
		join(BOWTIEDIR, "BAM_files/{sample}.bam"),
		join(BOWTIEDIR, "BAM_files/{sample}.bam.bai")
	params:
		outdir = join(BOWTIEDIR, "RM_files")
	output:
		join(BOWTIEDIR, "RM_files")
	shell:
		"if [ "$( wc -l {input} | cut -d' ' -f1)" -eq 0 ]; then; mv {input} {params.outdir}; fi"
'''

rule cut_contigs:
	input: 
		join(SPADESDIR, "contigs.fasta")
	output:
		join(CONCOCTDIR, "contigs_1k.bed"),
		join(CONCOCTDIR, "contigs_1k.fa")
	benchmark:
		join(BENCHDIR, "cutcontigs.benchmark.txt")
	shell:
		"cut_up_fasta.py {input} -c 1000 -o 0 --merge_last -b {output[0]} > {output[1]}"

rule cov_table:
	input:
		join(CONCOCTDIR, "contigs_1k.bed"),
		expand(join(BOWTIEDIR, "BAM_files/{sample}.bam.bai"), sample=IDS)
	output:
		join(CONCOCTDIR, "coverage_table.tsv")
	params:
		bam = expand(join(BOWTIEDIR, "BAM_files/{sample}.bam"), sample=IDS)
	benchmark:
		join(BENCHDIR, "covtable.benchmark.txt")
	shell:
		"concoct_coverage_table.py {input[0]} {params.bam} > {output}"

rule concoct:
	input:
		join(CONCOCTDIR, "contigs_1k.fa"),
		join(CONCOCTDIR, "coverage_table.tsv")
	output:
		join(CONCOCTDIR, "clustering_gt1000.csv")
	params:
		outdir = CONCOCTDIR
	benchmark:
		join(BENCHDIR, "concoct.benchmark.txt")
	shell:
		"~/.local/bin/concoct --threads 4 --composition_file {input[0]} --coverage_file {input[1]} -b {params.outdir}" 

rule merge_cutup:
	input:
		join(CONCOCTDIR, "clustering_gt1000.csv")
	output:
		join(CONCOCTDIR, "clustering_merged.csv")
	benchmark:
		join(BENCHDIR, "mergecutup.benchamrk.txt")
	shell:
		"merge_cutup_clustering.py {input} > {output}"

rule fasta_bins:
	input:
		join(CONCOCTDIR, "clustering_merged.csv")
	output:
		join(CONCOCTDIR, "fasta_bins/{bin}.fa")
	params:
		contigs = join(SPADESDIR, "contigs.fasta"),
		outdir = join(CONCOCTDIR, "fasta_bins")
	benchmark:
		join(BENCHDIR, "{bin}.fastabins.benchmark.txt")
	shell:
		"extract_fasta_bins.py {params.contigs} {input} --output_path {params.outdir}"

rule blast:
	input:
		join(CONCOCTDIR, "fasta_bins/{bin}.fa")
	output:
		join(BLASTDIR, "{bin}.out")
	params:
		blast_db = '~/anticrAss/new/blast/phage_genes'
	benchmark:
		join(BENCHDIR, "{bin}.blast.benchmark.txt")
	shell:
		"blastx -db {params.blast_db} -query {input} -out {output}"


rule readblast:
	input:
		join(BLASTDIR, "{bin}.out")
	output:
		join(BLASTDIR, "protein_hits/{bin}_blast.tsv")
	params:
		blast_reader = BLASTREADER,
		blastdir = BLASTDIR + "/"
	benchmark:
		join(BENCHDIR, "{bin}.readblast.benchmark.txt")
	shell:
		"set +o pipefail; mkdir blast/protein_hits/ | python {params.blast_reader} -i  {params.blastdir} -o {output}"




