# Work-in-progress for pipeline to take in SRA accession numbers and cross assemble for computational phage discovery. 

# STEP 1: Fastq-dump

#mkdir ~/anticrAss/new/run8/fastq

#while read line; do \
#echo "Downloading $line"; \
#/usr/local/genome/sra/current/bin/fastq-dump --outdir ~/anticrAss/new/run8/fastq -N 100001 -X 200000 --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line;\
#done <~/anticrAss/new/run8/SraAccList.txt

# STEP 1.5: Fastq Check (not working yet)

#while read line; do \
#if [ $(wc -l <~/anticrAss/new/run8/fastq/$line\_pass.fastq) -eq 400000 ]; then
#	echo "$line\_pass.fastq is complete";\
#else
#	echo "FAILED";\
#	rm ~/anticrAss/new/run8/fastq/$line\_pass.fastq
#fi
#done <~/anticrAss/new/run8/SraAccList.txt


# STEP 2: Combine into one large fastq file

echo "Concatenating files..."
cat ~/anticrAss/new/run8/fastq/*>~/anticrAss/new/run8/fastq/All_run8.fastq

# STEP 3: SPAdes Assembly

mkdir ~/anticrAss/new/run8/assembly
echo "Starting SPAdes Assembly..."
~/bin/anaconda3/envs/concoct_env/bin/python2 ~/anticrAss/old/SPAdes_Dir/SPAdes-3.5.0-Linux/bin/spades.py -m 1500 --12 ~/anticrAss/new/run8/fastq/All_run8.fastq --threads 16 -o ~/anticrAss/new/run8/assembly

# STEP 4: Create Bowtie2 index and scan the original reads against the assembled contigs

mkdir ~/anticrAss/new/run8/bowtie
mkdir ~/anticrAss/new/run8/bowtie/BAM_files
mkdir ~/anticrAss/new/run8/bowtie/SAM_files

echo "Creating Bowtie2 Index..."
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build ~/anticrAss/new/run8/assembly/contigs.fasta ~/anticrAss/new/run8/bowtie/All_run8idx



#####THIS CHUNK is the important parts for bowtie. Once the check is working, this is all we need
while read line; do \
	echo "Scanning $line with Bowtie2...";\
	~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q -x ~/anticrAss/new/run8/bowtie/All_run8idx -U ~/anticrAss/new/run8/fastq/$line\_pass*.fastq -S ~/anticrAss/new/run8/bowtie/SAM_files/$line.sam;\
	# STEP 4.5: SAM/BAM Check (works for correct files, still creates bam files of incorrect files. Try with only "if not" statement)

	if ~/bin/samtools-1.3.1/samtools quickcheck $line.sam ; then
		continue
	else
		echo "Converting sam to bam...";\
		~/bin/samtools-1.3.1/samtools view -bS ~/anticrAss/new/run8/bowtie/SAM_files/$line.sam | ~/bin/samtools-1.3.1/samtools sort -o ~/anticrAss/new/run8/bowtie/BAM_files/$line.bam;\
		~/bin/samtools-1.3.1/samtools index ~/anticrAss/new/run8/bowtie/BAM_files/$line.bam;\
	fi
done<~/anticrAss/new/run8/SraAccList.txt

# STEP 4.5 CORRECTION (can probably replace other step 4.5 but needs more testing; for now, okay to run both)
# This line will remove any empty bam files
for file in ~/anticrAss/new/run8/bowtie/BAM_files/*; do if [ "$( wc -l $file | cut -d' ' -f1)" -eq 0 ]; then rm "$file"; fi; done 


# STEP 5: CONCOCT Genome Binning

mkdir ~/anticrAss/new/run8/concoct_output

echo "Cutting up contigs for binning..."
~/bin/anaconda3/bin/python ~/bin/CONCOCT/scripts/cut_up_fasta.py ~/anticrAss/new/run8/assembly/contigs.fasta -c 1000 -o 0 --merge_last -b ~/anticrAss/new/run8/concoct_output/contigs_1k.bed > ~/anticrAss/new/run8/concoct_output/contigs_1k.fa

#echo "Creating coverage table..."
~/bin/anaconda3/bin/python ~/bin/CONCOCT/scripts/concoct_coverage_table.py ~/anticrAss/new/run8/concoct_output/contigs_1k.bed ~/anticrAss/new/run8/bowtie/BAM_files/*.bam > ~/anticrAss/new/run8/concoct_output/coverage_table.tsv

#echo "Running Concoct..."
~/.local/bin/concoct --threads 4 --composition_file ~/anticrAss/new/run8/concoct_output/contigs_1k.fa --coverage_file ~/anticrAss/new/run8/concoct_output/coverage_table.tsv -b ~/anticrAss/new/run8/concoct_output

python ~/bin/CONCOCT/scripts/merge_cutup_clustering.py ~/anticrAss/new/run8/concoct_output/clustering_gt1000.csv > ~/anticrAss/new/run8/concoct_output/clustering_merged.csv

mkdir ~/anticrAss/new/run8/concoct_output/fasta_bins

echo "Extracting genome bins..."
python  ~/bin/CONCOCT/scripts/extract_fasta_bins.py ~/anticrAss/new/run8/assembly/contigs.fasta ~/anticrAss/new/run8/concoct_output/clustering_merged.csv --output_path ~/anticrAss/new/run8/concoct_output/fasta_bins

mkdir ~/anticrAss/new/run8/blast/

for file in ~/anticrAss/new/run8/concoct_output/fasta_bins/*.fa; do
	echo $file
	blastx -db ~/anticrAss/new/blast/phage_genes -query $file -out ~/anticrAss/new/run8/blast/${file##*/}.00out;
	done

