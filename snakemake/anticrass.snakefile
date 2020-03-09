'''

Cross assembly pipeline for bacteriophage discovery from metagenomic samples.

To run this pipeline on antihll, use:

	snakemake -s anticrass.snakefile --cluster "qsub -cwd -V" --latency-wait 60

'''
from os.path import join


configfile: "anticrass.yaml"

IDS = []
with open("SraAccList.txt", 'r') as infile:
	for line in infile:
		IDS.append(line.strip())

print(IDS)

FASTQDIR = config["directories"]["fastq_files"]
SPADESDIR = config["directories"]["assembly"]
BOWTIEDIR = config["directories"]["bowtie"]
BLASTDIR = config["directories"]["blast_results"]

#FASTQ_IDX = range(1,3)

SPADES = config["commands"]["spades"]

rule all:
	input: 
		expand(join(FASTQDIR, "{sample}_pass_1.fastq"), sample=IDS), 
		join(FASTQDIR, "AllReads_1.fastq"),
		join(FASTQDIR, "AllReads_2.fastq"),
		join(SPADESDIR, "contigs.fasta"),
		expand(join(BOWTIEDIR, "contigs_idx.{index}.bt2"), index=range(1,5)),                                                                                                                  expand(join(BOWTIEDIR, "contigs_idx.rev.{index}.bt2"), index=range(1,3))

rule download:
	output:
		join(FASTQDIR, "{sample}_pass_1.fastq"),
		join(FASTQDIR, "{sample}_pass_2.fastq")
	params:
		sra_acc = "{sample}",
		outdir = FASTQDIR
	shell:
		"fastq-dump --outdir {params.outdir} -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {params.sra_acc}"


rule concat1:
	input:
		expand(join(FASTQDIR, "{sample}_pass_1.fastq"), sample=IDS)
	output:
		join(FASTQDIR, "AllReads_1.fastq")
	shell:
		"cat {input} > {output}"

rule concat2:
        input:
                expand(join(FASTQDIR, "{sample}_pass_2.fastq"), sample=IDS)
        output:
                join(FASTQDIR, "AllReads_2.fastq")
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
	shell:
		"bowtie2-build {input} {params.idx}" 

