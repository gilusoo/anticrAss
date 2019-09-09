# This is a pipeline that takes in a list of SRR's in a txt file and returns predicted genome bins
# The commands should work for all runs within the sra bioproject PRJNA327423

# STEP 1: Download the necessary datasets from the SRA using fastq-dump
# (This will likely be replaced with fasterq-dump in the future.)

mkdir ~/anticrAss/new/run3/fastq

while read line; do \
echo "Downloading $line"; \
fastq-dump --outdir ~/anticrAss/new/run3/fastq --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line;\
done <~/anticrAss/new/run3/SraAccList.txt

# STEP 2: Concat all fastq files into one for big assembly

cat ~/anticrAss/new/run3/fastq/*>~/anticrAss/new/run3/fastq/All_run3.fastq

# STEP 3: Run one big spades assembly on the fastq file we created

mkdir ~/anticrAss/new/run3/assembly
~/anticrAss/old/SPAdes_Dir/SPAdes-3.5.0-Linux/bin/spades.py -m 1000 --12 ~/anticrAss/new/run3/fastq/All_run3.fastq -o ~/anticrAss/new/run3/assembly

# STEP 4: Create Bowtie2 index and scan the original files against the assembled contigs

mkdir ~/anticrAss/new/run3/bowtie
echo "Creating Bowtie2 Index..."
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build ~/anticrAss/new/run3/assembly/contigs.fasta ~/anticrAss/new/run3/bowtie/All_run3idx

mkdir ~/anticrAss/new/run3/bowtie/BAM_files
mkdir ~/anticrAss/new/run3/bowtie/SAM_files

while read line; do \
echo "Scanning $line with Bowtie2...";\
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q -x ~/anticrAss/new/run3/bowtie/All_run3idx -U ~/anticrAss/new/run3/fastq/$line\_pass.fastq -S ~/anticrAss/new/run3/bowtie/SAM_files/$line.sam;\
echo "Converting sam to bam...";\
~/bin/samtools-1.3.1/samtools view -bS ~/anticrAss/new/run3/bowtie/SAM_files/$line.sam | ~/bin/samtools-1.3.1/samtools sort -o ~/anticrAss/new/run3/bowtie/BAM_files/$line.sorted.bam;\
~/bin/samtools-1.3.1/samtools index ~/anticrAss/new/run3/bowtie/BAM_files/$line.sorted.bam;\
done<~/anticrAss/new/run3/SraAccList.txt

echo "Scanning with big fastq..."
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q -x ~/anticrAss/new/run3/bowtie/All_run3idx -U ~/anticrAss/new/run3/fastq/All_run3.fastq -S ~/anticrAss/new/run3/bowtie/SAM_files/All_run3.sam
echo "Converting sam to bam..."
~/bin/samtools-1.3.1/samtools views -bS ~/anticrAss/new/run3/bowtie/SAM_files/All_run3.sam | ~/bin/samtools-1.3.1/samtools sort -o ~/anticrAss/new/run3/bowtie/BAM_files/All_run3.sorted.bam
~/bin/samtools-1.3.1/samtools index ~/anticrAss/new/run3/bowtie/BAM_files/All_run3.sorted.bam


# STEP 5: CONCOCT Genome Binning

mkdir ~/anticrAss/new/run3/concoct_output

echo "Cutting up contigs for binning..."
ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/cut_up_fasta.py ~/anticrAss/new/run3/assembly/contigs.fasta -c 1000 -o 0 --merge_last -b ~/anticrAss/new/run3/concoct_output/contigs_1k.bed > ~/anticrAss/new/run3/concoct_output/contigs_1k.fa

echo "Creating coverage table..."
ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/concoct_coverage_table.py ~/anticrAss/new/run3/concoct_output/contigs_1k.bed ~/anticrAss/new/run3/bowtie/BAM_files/All_run3.sorted.bam > ~/anticrAss/new/run3/concoct_output/coverage_table.tsv

echo "Running Concoct..."
ssh tatabox ~/bin/anaconda3/envs/concoct_env/bin/concoct --composition_file ~/anticrAss/new/run3/concoct_output/contigs_1k.fa --coverage_file ~/anticrAss/new/run3/concoct_output/coverage_table.tsv -b ~/anticrAss/new/run3/concoct_output

ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/merge_cutup_clustering.py ~/anticrAss/new/run3/concoct_output/clustering_gt1000.csv > ~/anticrAss/new/run3/concoct_output/clustering_merged.csv

mkdir ~/anticrAss/new/run3/concoct_output/fasta_bins

echo "Extracting genome bins..."
ssh tatabox python  ~/bin/anaconda3/envs/concoct_env/bin/extract_fasta_bins.py ~/anticrAss/new/run3/assembly/contigs.fasta ~/anticrAss/new/run3/concoct_output/clustering_merged.csv --output_path ~/anticrAss/new/run3/concoct_output/fasta_bins
