import re
import argparse
import sys
import os

class Blast_Reader:
    """
    This is a class containing functions to help parse through results files created by NCBI's BLAST for the command line.
    """
    def __init__(self, path):
        """
        Initializes the Blast_Reader class.
        :param path: The path to the directory of output file from BLAST on the command line.
        """
        self.name = path

    def hit_summary(self) -> list:
        """
        Creates a list of lists, where each inner list contains the information provided in the hit summary.
        :return: A list of lists
        """
        summary = [line.strip().split() for line in open(self.name) if re.match(r"[A-Z]{3}[0-9]{4}", line[0:8])]
        return summary

    def hit_cov_e(self, summary) -> dict:
        """
        Creates a dictionary in which each key is a protein accession number and each value is a list containing the
        percent coverage of that hit, as well as e-value.
        :param summary: The list of lists created by the hit_summary method.
        :return: A dictionary
        """
        prot_to_cov_e = {}
        for line in summary:
            prot_to_cov_e[line[0]] = [line[-2], line[-1]]
        return prot_to_cov_e

    def write_hits(self, e_dict: dict, bin, outfile) -> None:
        with open(outfile, 'w') as outfile:
            outfile.write('Bin,Protein,% Coverage,e-Value\n')
            lines = [f'{bin},{k},{v[0]},{v[1]}\n' for k, v in e_dict.items() if float(v[0]) >= 40.0 and
                                                      (len(re.findall(r"e-0[3456789]", v[1])) >= 1 or
                                                      len(re.findall(r"e-[123456789]{2}", v[1])) >= 1)]
            for line in lines:
                outfile.write(line)
        return

    # ----------------------------------------------------

    def get_contig(self, line):
        """
        Pulls the name of the current contig being analyzed.
        :param line: A line of the BLAST results file
        :return: string, contig name
        """
        return line.split()[1]

    def prot_dict(self, contig, line):
        """
        Creates a dictionary with hits of scores greater than 40% coverage
        :param line: A line of the BLAST results file
        :return: dictionary in the format: {contig: {protein: (score, e-value) ...} ...}
        """
        prot_dict = {}
        if float(line.split()[-2]) >= 40.0:
            prot_dict[contig] = {line.split()[0]: (line.split()[-2], line.split()[-1])}
        return prot_dict

    def write_prots(self, contig, prot_dict, outfile):
        """
        Writes the resulting protein hits to a csv file
        :param contig: The name of the contig being analyzed
        :param prot_dict: A dictionary of contigs to proteins to score, e-value
        :param outfile: The file to write to
        :return:
        """
        with open(outfile, 'a') as outfile:
            for k, v in prot_dict.items():
                for prot, nums in v.items():
                    outfile.write(f'{contig}\t{prot}\t{nums[0]}\t{nums[1]}\n')
        return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Directory of BLAST results files')
    parser.add_argument('-o', '--output', help='Output csv file to write to')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    br = Blast_Reader(args.input)

    for root, dirs, files in os.walk(args.input):
        for name in files:
            if root == args.input:
                file = os.path.join(root, name)
                output = root + f'protein_hits/' + str(name.split('.')[0]) + '_blast.tsv'
                with open(output, 'w') as outfile:
                    outfile.write('Contig\tProtein\tScore\te-Value\n')
                for line in open(file):
                    if line.startswith('Query= '):
                        curr_contig = br.get_contig(line)
                    if re.match(r"[A-Z]{3}[0-9]{4}", line[0:8]):
                        #if float(line.split()[-2]) >= 40.0:
                        if float(line.split()[-1]) <= 0.50:
                            br.write_prots(curr_contig, br.prot_dict(curr_contig, line), output)
