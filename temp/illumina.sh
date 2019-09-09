# Work-in-progress for pipeline to take in SRA accession numbers and cross assemble for computational phage discovery. 

# STEP 1: Fastq-dump

while read line; do \
echo "Downloading $line"; \
fastq-dump --outdir ~/anticrAss/new/run4/fastq --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line;\
done <~/anticrAss/new/run4/SraAccList.txt

# STEP 2: Combine into one large fastq file

cat ~/anticrAss/new/run4/fastq/*>~/anticrAss/new/run4/fastq/All_run4.fastq

# STEP 3: SPAdes Assembly

mkdir ~/anticrAss/new/run4/assembly
python2 ~/anticrAss/old/SPAdes_Dir/SPAdes-3.5.0-Linux/bin/spades.py -m 1000 --12 ~/anticrAss/new/run4/fastq/All_run4.fastq -o ~/anticrAss/new/run4/assembly

# STEP 4: Create Bowtie2 index and scan the original reads against the assembled contigs

echo "Creating Bowtie2 Index..."
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build ~/anticrAss/new/run4/assembly/contigs.fasta ~/anticrAss/new/run4/bowtie/All_run4idx

while read line; do \
echo "Scanning $line with Bowtie2...";\
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q -x ~/anticrAss/new/run4/bowtie/All_run4idx -U ~/anticrAss/new/run4/fastq/$line\_pass.fastq -S ~/anticrAss/new/run4/bowtie/SAM_files/$line.sam;\
echo "Converting sam to bam...";\
~/bin/samtools-1.3.1/samtools view -bS ~/anticrAss/new/run4/bowtie/SAM_files/$line.sam | ~/bin/samtools-1.3.1/samtools sort -o ~/anticrAss/new/run4/bowtie/BAM_files/$line.sorted.bam;\
~/bin/samtools-1.3.1/samtools index ~/anticrAss/new/run4/bowtie/BAM_files/$line.sorted.bam;\
done<~/anticrAss/new/run4/SraAccList.txt

echo "Scanning with big fastq..."
~/bin/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q -x ~/anticrAss/new/run4/bowtie/All_run4idx -U ~/anticrAss/new/run4/fastq/All_run4.fastq -S ~/anticrAss/new/run4/bowtie/SAM_files/All_run4.sam
echo "Converting sam to bam..."
~/bin/samtools-1.3.1/samtools views -bS ~/anticrAss/new/run4/bowtie/SAM_files/All_run4.sam | ~/bin/samtools-1.3.1/samtools sort -o ~/anticrAss/new/run4/bowtie/BAM_files/All_run4.sorted.bam
~/bin/samtools-1.3.1/samtools index ~/anticrAss/new/run4/bowtie/BAM_files/All_run4.sorted.bam

# STEP 5: CONCOCT Genome Binning

echo "Cutting up contigs for binning..."
ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/cut_up_fasta.py ~/anticrAss/new/run4/assembly/contigs.fasta -c 1000 -o 0 --merge_last -b ~/anticrAss/new/run4/concoct_output/contigs_1k.bed > ~/anticrAss/new/run4/concoct_output/contigs_1k.fa

echo "Creating coverage table..."
ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/concoct_coverage_table.py ~/anticrAss/new/run4/concoct_output/contigs_1k.bed ~/anticrAss/new/run4/bowtie/BAM_files/*.sorted.bam > ~/anticrAss/new/run4/concoct_output/coverage_table.tsv

echo "Running Concoct..."
ssh tatabox ~/bin/anaconda3/envs/concoct_env/bin/concoct --composition_file ~/anticrAss/new/run4/concoct_output/contigs_1k.fa --coverage_file ~/anticrAss/new/run4/concoct_output/coverage_table.tsv -b ~/anticrAss/new/run4/concoct_output

ssh tatabox python ~/bin/anaconda3/envs/concoct_env/bin/merge_cutup_clustering.py ~/anticrAss/new/run4/concoct_output/clustering_gt1000.csv > ~/anticrAss/new/run4/concoct_output/clustering_merged.csv

mkdir ~/anticrAss/new/run4/concoct_output/fasta_bins

echo "Extracting genome bins..."
ssh tatabox python  ~/bin/anaconda3/envs/concoct_env/bin/extract_fasta_bins.py ~/anticrAss/new/run4/assembly/contigs.fasta ~/anticrAss/new/run4/concoct_output/clustering_merged.csv --output_path ~/anticrAss/new/run4/concoct_output/fasta_bins
