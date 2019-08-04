"""
Reads in a text file containing a list of SRA accessions. Then goes through one run at a time
to collect the available metadata and writes the desired fields to a new file.
"""

import sys
import requests
from bs4 import BeautifulSoup
import argparse


def sra_acc(infile):
    acc_list = []
    with open(infile, 'r'):
        for line in infile:
            line = line.strip()
            acc_list.append(line)
    return acc_list

def pull_metadata(list, outfile, *args):
    header = ['Run', args]
    outlines = []
    for run in list:
        req = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={}'.format(run))
        soup = BeautifulSoup(req.text, 'lxml')
        line = soup.find("div", {'class': 'ph experiment'})
        all_attributes = line.find_all('th')
        all_att_values = line.find_all('td')
        attributes_dict = {}
        for idx, att in enumerate(all_attributes):
            attributes_dict[att] = all_att_values[idx].contents[0]
        currline = 'Run'
        for item in args:
            data = attributes_dict.get(item, 'N/A')
            currline = currline + ',' + data
        outlines.append(currline)
        currline = 'Run'
    with open(outfile, 'w'):
        for item in header:
            outfile.write('{},'.format(item))
        outfile.write('\n')
        for line in currline:
            outfile.write('{}\n'.format(line))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input SRA accession list')
    parser.add_argument('-o', '--output', help='Output csv file to write to')
    parser.add_argument('-a', '--attributes', nargs='+', help='List of desired attributes from each run; '
                                                              'separated by a single space')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)


acclist = sra_acc(args.input)
pull_metadata(acclist, args.output, args.attributes)
