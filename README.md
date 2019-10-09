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
