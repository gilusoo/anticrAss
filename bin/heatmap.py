import os
import sys
import pandas as pd
import argparse
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns


class Heatmaps:
    '''
    This is a class containing functions necessary for the generation of a heatmap from read mapping data. The bam files
    and indexed bam files should first be run through sam_hits.py to generate a "sam_hits" or "bam_hits" file to use as
    input.
    # Later can be edited to generate a heatmap from bam/bam.bai files in a directory (see Sam_Reader init)
    '''
    def __init__(self, path):
        self.name = path

    def get_contigs(self, fasta):
        '''
        Extracts the names of the contigs from the headers of a fasta file and makes them into a list.
        :return: A dictionary of contig names to their index number. A list of just the contig names can be extracted
        from this dictionary.
        '''
        contig_idx = {}
        with open(fasta) as infile:
            contigs = [x.strip().replace('>','') for x in infile if x.startswith('>')]
            for idx, contig in enumerate(contigs):
                contig_idx[str(contig)] = idx
        return contig_idx

    def hits_dict(self, data, contigs):
        '''
        Creates a dictionary of the hit data in the format {run: {contig: hits, ...}, ...}
        :param data: Input file of bam_hits data (can be generated using sam_hits.py)
        :param contigs: A list of contig names to be analyzed
        :return: Dictionary
        '''
        runs_to_hits = {}
        with open(data) as infile:
            next(infile)
            for line in infile:
                line = line.strip().split('\t')
                num_hits = int(line[3])
                if num_hits == 0:
                    continue
                run = os.path.basename(line[0]).split('.')[0]
                runs_to_hits[run] = runs_to_hits.get(run, {contig: 0  for contig in contigs})
                hit = line[1]
                runs_to_hits[run][hit] = num_hits
        return runs_to_hits

    def group_runs(self, annot):
        '''
        Creates a dictionary grouping the SRA run accessions based on their source.
        :param annot: Raw SRA run annotations file
        :return: Dictionary
        '''
        cat_to_runs = {}
        with open(annot, encoding='utf-8', errors='ignore') as infile:
            next(infile)
            for line in infile:
                line = line.strip().split('\t')
                if len(line) > 4:
                    runs = line[5].strip().split(',')
                if len(line) > 3 and line[3] != '':
                    cat = line[3]
                    cat_to_runs[cat] = cat_to_runs.get(cat, []) + runs
        return cat_to_runs

    def heatmap_df(self, hits_dict, contigs, cat_dict, groups):
        '''
        Creates a pandas df in the correct format for a heatmap.
        :param hits_dict: Dictionary in the format of {run: {contig: hits, ...}, ...}
        :param contigs: List of contig names to be analyzed
        :param cat: List of runs with a specific annotation
        :return: Pandas dataframe ready to create heatmap
        '''
        df = []
        yaxis = []
        lengths = []
        count = 0
        for cat in groups:
            runs = cat_dict[cat]
            cat_count = 0

            def custom_sort(k):
                """
                this definition takes in one inner element of whatever is being sorted
                based on this inner element, this def should return an int
                """
                innerd = hits_dict[k]
                return sum(innerd.values())

            for run in sorted(hits_dict.keys(), key=custom_sort, reverse=True):
                # print(run, sum(hits_dict[run].values()))
                if run in runs:
                    if cat_count < 1000:
                        count += 1
                        cat_count =+ 1
                        yaxis.append(run)

                    # print(run, sum(hits_dict[run].values()))
                        curr = []
                        for contig in contigs[:20]:
                            # if hits_dict[run][contig] < 150:
                            curr.append(hits_dict[run][contig])
                            # else:
                            #     curr.append(150)
                        df.append(curr)
            lengths.append(count)
        for x in yaxis[:500]:
            print(x)
        return pd.DataFrame(df, index=yaxis, columns=contigs), lengths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--fasta', help='Input fasta file of genome')
    parser.add_argument('-d', '--data', help='Input tsv file containing mapping data. (Obtained using sam_hits.py)')
    parser.add_argument('-o', '--output', help='Output filename and path for the resulting heatmap')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    hm = Heatmaps(args.fasta)

    test = hm.get_contigs(args.fasta)
    # print(test)
    # print(len(test))
    # print(hm.hits_dict(args.data, hm.get_contigs(args.fasta)))
    # print(hm.heatmap_df(hm.hits_dict(args.data, hm.get_contigs(args.fasta)), hm.get_contigs(args.fasta)))

    # print(hm.group_runs('../sra_annotations.tsv'))      # The annotations file must be in same directory as the script

    ##############

    groups = hm.group_runs(r'sra_annotations.tsv')
    of_interest = ['human gut', 'human respiratory', 'human skin', 'animal', 'marine deep', 'marine surface', 'soil']

    data = hm.heatmap_df(hm.hits_dict(args.data, hm.get_contigs(args.fasta)), hm.get_contigs(args.fasta),
                         hm.group_runs('sra_annotations.tsv'), of_interest)[0]
    breaks = hm.heatmap_df(hm.hits_dict(args.data, hm.get_contigs(args.fasta)), hm.get_contigs(args.fasta),
                         hm.group_runs('sra_annotations.tsv'), of_interest)[1]

    heat = sns.heatmap(data, cmap='inferno', vmin=0, vmax=150)
    plt.show(block=True)
    plt.hlines(breaks, xmin=0, xmax=len(hm.get_contigs(args.fasta)), colors='w')
    plt.savefig(args.output)

