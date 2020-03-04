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
BLASTDIR = config["directories"]["blast_results"]

rule all:
	input: 
		expand(join(FASTQDIR, "{sample}_pass.fastq"), sample=IDS),
		join(FASTQDIR, "AllReads.fastq")

rule download:
	output:
		join(FASTQDIR, "{sample}_pass.fastq")
	params:
		sra_acc = "{sample}"
		#outdir = os.path(FASTQDIR)
	shell:
		"fastq-dump --outdir fastq/ -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {params.sra_acc}"
		#fastq-dump --outdir ~/anticrAss/new/run8/fastq -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line

rule concat:
	input:
		expand(join(FASTQDIR, "{sample}_pass.fastq"), sample=IDS)
	output:
		join(FASTQDIR, "AllReads.fastq")
	params:
		fastq = "os.path(FASTQDIR)"
	shell:
		"cat {input} > {output}"

'''
rule assemble:
	output:
	params:
	shell:
		"python2 spades.py -m 1500 --12
		~/bin/anaconda3/envs/concoct_env/bin/python2 ~/anticrAss/old/SPAdes_Dir/SPAdes-3.5.0-Linux/bin/spades.py -m 1500 --12 ~/anticrAss/new/run8/fastq/All_run8.fastq --threads 16 -o ~/anticrAss/new/run8/assembly
'''

#BINS, = glob_wildcards(join(BLASTDIR, "{bin}.out"))





