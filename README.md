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
        
   After all of the datasets have been successfully downloaded, we need to concatenate them into a single FASTA file for assembly. This can be done using the following command:
   
     cat /anticrAss/experiment/fastq/* > /anticrAss/experiment/fastq/All_experiment.fastq
        
## De Novo Assembly of Raw Sequence Reads

  Assembling raw sequence reads is essential for making use of the entirity of information the reads contain. For novel sequence discovery, reference-guided assembly techniques are insufficient. Instead, we use SPAdes, an assembly pipeline that works independent of a reference genome by assembling contigs based on overlapping reads. 
  
   The following is a good starting point for a SPAdes assembly command. Documentation of all command line options can be found at their [github repo](https://github.com/ablab/spades).
   
     mkdir ~/anticrAss/new/run4/assembly
     python SPAdes-3.5.0-Linux/bin/spades.py -m 1000 --12 /anticrAss/experiment/fastq/All_run4.fastq -o /anticrAss/experiment/assembly
     
