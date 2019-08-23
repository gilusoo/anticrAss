# This Makefile contains multiple variations of the anticrAss pipeline to be used for downloading and analyzing
# data from the SRA.
# The commands can be run from the terminal by typing "make command_name". EX. make setup

# ---------------SETUP---------------

setup:
	mkdir -p Experiment/assembly Experiment/fastq Experiment/concoct_output
	mkdir -p Experiment/bowtie/BAM_files Experiment/bowtie/SAM_files

# ---------DATA SELECTION------------

illumina_gut_single:
	python gen_acc_list.py -i data/gut/gut_db.txt -o Experiment -n 10 -p Illumina -l Single

illumina_gut_paired:
	python gen_acc_list.py -i data/gut/gut_db.txt -o Experiment -n 10 -p Illumina -l Paired
