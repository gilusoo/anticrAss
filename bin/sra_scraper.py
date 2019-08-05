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
    with open(infile, 'r') as infile:
        for line in infile:
            line = line.strip()
            acc_list.append(line)
    return acc_list

def pull_metadata(list, outfile, att_list):
    header = 'Run,{}'.format(','.join(att_list))
    outlines = []
    for run in list:
        req = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={}'.format(run))
        soup = BeautifulSoup(req.text, 'lxml')
        line = soup.find("div", {'class': 'ph experiment'})
        all_attributes = line.find_all('th')
        all_att_values = line.find_all('td')
        attributes_dict = {}
        for idx, att in enumerate(all_attributes):
            attributes_dict[att.contents[0]] = all_att_values[idx].contents[0]
        currline = run
        for item in att_list:
            data = attributes_dict.get(item, 'N/A')
            currline = currline + ',' + data
        outlines.append(currline)
    with open(outfile, 'w') as outfile:
        outfile.write(header)
        outfile.write('\n')
        for line in outlines:
            outfile.write('{}\n'.format(line))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input SRA accession list')
    parser.add_argument('-o', '--output', help='Output csv file to write to')
    parser.add_argument('-a', '--attributes', nargs='+', help='List of desired attributes from each run;'
                                                              'format: \'[<att1>,<att2>,...])\'')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    acclist = sra_acc(args.input)
    pull_metadata(acclist, args.output, args.attributes)
