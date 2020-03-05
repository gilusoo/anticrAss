'''

Cross assembly pipeline for bacteriophage discovery from metagenomic samples.

To run this pipeline on antihll, use:

	snakemake -s anticrass.snakefile --cluster "qsub -cwd -V" --latency-wait 60

'''
from os.path import join


configfile: "anticrass.yaml"

IDS = []
with open("SraAccList.txt", 'r') as infile:
	for acc in infile:
		IDS.append(acc.strip())

FASTQDIR = config["directories"]["fastq_files"]
SPADESDIR = config["directories"]['assembly']
BLASTDIR = config["directories"]["blast_results"]

SPADES = config["commands"]["spades"]

rule all:
	input: 
		expand(join(FASTQDIR, "{sample}_pass.fastq"), sample=IDS),
		join(FASTQDIR, "AllReads.fastq"),
		join(SPADESDIR, "contigs.fasta")

rule download:
	output:
		join(FASTQDIR, "{sample}_pass.fastq")
	params:
		sra_acc = "{sample}",
		outdir = "os.path(FASTQDIR)"
	shell:
		"fastq-dump --outdir {params.outdir} -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {params.sra_acc}"
		#fastq-dump --outdir ~/anticrAss/new/run8/fastq -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line

rule concat:
	input:
		expand(join(FASTQDIR, "{sample}_pass.fastq"), sample=IDS)
	output:
		join(FASTQDIR, "AllReads.fastq")
	shell:
		"cat {input} > {output}"

rule assemble:
	input:
		join(FASTQDIR, "AllReads.fastq")
	output:
		join(SPADESDIR, "contigs.fasta"),
		join(SPADESDIR, "contigs.fastg"),
		join(SPADESDIR, "dataset.info"),
		join(SPADESDIR, "scaffolds.fasta"),
		join(SPADESDIR, "scaffolds.fastg"), 
		join(SPADESDIR, "spades.log")
	params:
		spades_assemble = SPADES,
		outdir = config["directories"]["assembly"]
	shell:
		"python2 {params.spades_assemble} -m 1500 --12 {input} --threads 16 --only-assembler -o {params.outdir}"
		#~/bin/anaconda3/envs/concoct_env/bin/python2 ~/anticrAss/old/SPAdes_Dir/SPAdes-3.5.0-Linux/bin/spades.py -m 1500 --12 ~/anticrAss/new/run8/fastq/All_run8.fastq --threads 16 -o ~/anticrAss/new/run8/assembly


#BINS, = glob_wildcards(join(BLASTDIR, "{bin}.out"))





