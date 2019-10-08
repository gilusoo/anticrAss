echo "Cutting up contigs for binning..."
python cut_up_fasta.py [PATH TO contigs.fasta] -c 1000 -o 0 --merge_last -b [PATH TO BED FILE /concoct_output/contigs_1k.bed] > [PATH TO OUTPUT FILE /concoct_output/contigs_1k.fa]

echo "Creating coverage table..."
python concoct_coverage_table.py [PATH TO BED FILE /concoct_output/contigs_1k.bed] [PATH TO BAM FILES /BAM_files/*.sorted.bam] > [PATH TO OUTPUT FILE /concoct_output/coverage_table.tsv]

echo "Running Concoct..."
concoct --composition_file [PATH TO CUT UP CONTIGS /concoct_output/contigs_1k.fa] --coverage_file [PATH TO COVERAGE TABLE /concoct_output/coverage_table.tsv] -b [PATH TO OUTPUT DIR /concoct_output]

python merge_cutup_clustering.py [PATH TO CSV FILE /concoct_output/clustering_gt1000.csv] > [PATH TO OUTPUT FILE /concoct_output/clustering_merged.csv]

mkdir [PATH TO DIR /concoct_output/fasta_bins]

echo "Extracting genome bins..."
python extract_fasta_bins.py [PATH TO ASSEMBLED CONTIGS /assembly/contigs.fasta] [PATH TO MERGED CSV /concoct_output/clustering_merged.csv] --output_path [PATH TO BINS DIR /concoct_output/fasta_bins]
