# anticrAss: A Computational Pipeline for Novel Phage Discovery

The 2014 computational discovery of crAssphage in the human gut sparked the idea that several other undocumented bacteriophages may be living in the same environment undetected. Here, we introduce anticrAss, a pipeline comprised of computational tools to analyze the species composition of metagenomes independent of a reference database. 

## Data Selection

The first step to making use of this pipeline is to create an SRA Accession List of metagenomes to be analyzed. There are a few methods for this. 

1. Create your own SRA Accession list using the SRA Run Selector. 
  
    Simply check the box next to any sample that is of interest. When you have finsihed selecting all of the datasets that you plan to use, download the accession list for you selected runs as SraAccList.txt
  
   **It is important that your accession list contains runs of the same format. You may consider using the SRA Run Selector filter options for:      
      - Platform
      - Strategy
      - LibraryLayout

2. Create a random SRA Accession List from one of the provided databases of metagenomes. 

3. Generate an SRA Accession List of runs with consistent data using web scraping techniques.

    A python script has been provided for creating an SRA Accession List of the desired platform (Illumina, Ion-torrent), library layout (Single, Paired), and sample source (Soil, Fresh Water, Marine, or Gut). This is currently the only method for creating an accession list from the gut database provided. All other databases can randomly generate an accession list without requiring a web scrape. 
    
    Here is an example of the command needed to generate an accesssion list using this method:
    
        python3 gen_acc_list.py -i /data/gut/gut_db.txt -n 10 -p Illumina -l Single -o /anticrAss/experiment/SraAccList.txt
    
    The above command would create an accession list of 10 gut metagenomes that are Illumina single reads. 
        
## Downloading SRA Datasets

  Once you have a text file containing the SRA accessions of all runs of interest, you are ready to obtain the actual files. This can be done fairly easily using fastq-dump. 
  
   The downloads can be done in a loop using the following code:
   
    while read line; do \
        echo "Downloading $line"; \
        fastq-dump --outdir /anticrAss/experiment/fastq --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line;\
        done <~/anticrAss/experiment/SraAccList.txt
        
   After all of the datasets have been successfully downloaded, we need to concatenate them into a single FASTQ file for assembly. This can be done using the following command:
   
     cat /anticrAss/experiment/fastq/* > /anticrAss/experiment/fastq/All_experiment.fastq
        
## De Novo Assembly of Raw Sequence Reads

  Assembling raw sequence reads is essential for making use of the entirity of information the reads contain. For novel sequence discovery, reference-guided assembly techniques are insufficient. Instead, we use SPAdes, an assembly pipeline that works independent of a reference genome by assembling contigs based on overlapping reads. 
  
   The following is a good starting point for a SPAdes assembly command. Documentation of all command line options can be found at their [github repo](https://github.com/ablab/spades).
   
     mkdir ~/anticrAss/new/run4/assembly
     python SPAdes-3.5.0-Linux/bin/spades.py -m 1000 --12 /anticrAss/experiment/fastq/All_experiment.fastq -o /anticrAss/experiment/assembly
     
## Read Mapping to Assembled Contigs

  Assembled contigs must be verified by mapping the raw reads to the assembly to check for coverage. The information contained in the results of this read mapping will also be used for genome binning. This is all done using Bowtie2.
  
  We must first create an index of the assembled contigs for the reads to be mapped against. This can be done using the following command:
  
      echo "Creating Bowtie2 Index..."
      bowtie2-build /anticrAss/experiment/assembly/contigs.fasta /anticrAss/experiment/bowtie/All_experimentidx
      
  The index is now ready for read mapping. The results of this will be in the form of a sam file. We will want to sort and convert these files from sam (human-readable) to bam (binary) for future use. This can all be done in one loop using the following commands for bowtie2 and samtools:
  
    while read line; do \
        echo "Scanning $line with Bowtie2...";\
        bowtie2 -q -x /anticrAss/experiment/bowtie/All_experimentidx -U /anticrAss/experiment/fastq/$line\_pass.fastq -S /anticrAss/experiment/bowtie/SAM_files/$line.sam;\
        echo "Converting sam to bam...";\
        samtools view -bS /anticrAss/experiment/bowtie/SAM_files/$line.sam | samtools sort -o ~/anticrAss/experiment/bowtie/BAM_files/$line.sorted.bam;\
        samtools index /anticrAss/experiment/bowtie/BAM_files/$line.sorted.bam;\
        done<~/anticrAss/experiment/SraAccList.txt
        
  Finally, we must also create sorted bam files for our large FASTQ file for future use. This can be done using the following command:
  
      echo "Scanning with big fastq..."
      bowtie2 -q -x /anticrAss/experiment/bowtie/All_experimentidx -U /anticrAss/experiment/fastq/All_experiment.fastq -S /anticrAss/experiment/bowtie/SAM_files/All_experiment.sam
      echo "Converting sam to bam..."
      samtools views -bS /anticrAss/experiment/bowtie/SAM_files/All_experiment.sam | samtools sort -o /anticrAss/experiment/bowtie/BAM_files/All_experiment.sorted.bam
      samtools index /anticrAss/experiment/bowtie/BAM_files/All_experiment.sorted.bam
      
## Genome Binning

  Now that we have a collection of assembled contigs and read mapping information in bam files, we can attempt to sort the contigs into predicted genome bins. This is done using CONCOCT, a binning software that sorts contigs using two methods:
  
  - Coverage
  - Composition
      
  The general method by which this program is able to accomplish this task is outlined in the following steps.
  
  1. Cutting Up Contigs
      
          python cut_up_fasta.py /anticrAss/experiment/assembly/contigs.fasta -c 1000 -o 0 --merge_last -b /anticrAss/experiment/concoct_output/contigs_1k.bed > /anticrAss/experiment/concoct_output/contigs_1k.fa
          
  2. Creating a Coverage Table
  
          python concoct_coverage_table.py /anticrAss/experiment/concoct_output/contigs_1k.bed /anticrAss/experiment/bowtie/BAM_files/*.sorted.bam > /anticrAss/experiment/concoct_output/coverage_table.tsv
          
  3. Running CONCOCT Binning
  
          concoct --composition_file ~/anticrAss/experiment/concoct_output/contigs_1k.fa --coverage_file /anticrAss/experiment/concoct_output/coverage_table.tsv -b /anticrAss/experiment/concoct_output
          
  4. Merging Cut Up Clusters
  
          python merge_cutup_clustering.py /anticrAss/experiment/concoct_output/clustering_gt1000.csv > /anticrAss/experiment/concoct_output/clustering_merged.csv
          
  5. Extracting Genome Bins
  
          mkdir ~/anticrAss/experiment/concoct_output/fasta_bins
          python extract_fasta_bins.py /anticrAss/experiment/assembly/contigs.fasta /anticrAss/experiment/concoct_output/clustering_merged.csv --output_path /anticrAss/experiment/concoct_output/fasta_bins
