# This Makefile contains multiple variations of the anticrAss pipeline to be used for downloading and analyzing
# data from the SRA.
# The commands can be run from the terminal by typing "make command_name". EX. make setup

# ---------------SETUP---------------

setup:
	mkdir -p Experiment/assembly Experiment/fastq Experiment/concoct_output
	mkdir -p Experiment/bowtie/BAM_files Experiment/bowtie/SAM_files

# ---------DATA SELECTION------------

illumina_gut_single:
	python bin/gen_acc_list.py -i data/gut/gut_db.txt -o Experiment -n 10 -p Illumina -l Single

illumina_gut_paired:
	python bin/gen_acc_list.py -i data/gut/gut_db.txt -o Experiment -n 10 -p Illumina -l Paired

illumina_soil_single:
	python bin/gen_acc_list.py -i data/soil/soil_db_illumina.txt -o Experiment -n 10 -p Illumina -l Single

illumina_soil_paired:
	python bin/gen_acc_list.py -i data/soil/soil_db_illumina.txt -o Experiment -n 10 -p Illumina -l Paired

illumina_marine_single:
	python bin/gen_acc_list.py -i data/marine/marine_db_illumina.txt -o Experiment -n 10 -p Illumina -l Single

illumina_marine_paired:
	python bin/gen_acc_list.py -i data/marine/marine_db_illumina.txt -o Experiment -n 10 -p Illumina -l Paired

illumina_fresh_water_single:
	python bin/gen_acc_list.py -i data/fresh_water/fresh_water_db_illumina.txt -o Experiment -n 10 -p Illumina -l Single

illumina_fresh_water_paired:
	python bin/gen_acc_list.py -i data/fresh_water/fresh_water_db_illumina.txt -o Experiment -n 10 -p Illumina -l Paired
	
