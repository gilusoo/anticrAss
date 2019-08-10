
import sys
import requests
from bs4 import BeautifulSoup
import argparse

def sra_acc(infile):
    with open(infile, 'r') as infile:
       acc_list = [line.strip() for line in infile]
    return acc_list

def pull_metadata(acc_list, attributes, outfile):
    outlines = []
    outfile = outfile + '_meta.csv'
    for run in acc_list:
        req = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={}'.format(run))
        soup = BeautifulSoup(req.text, 'lxml')
        line = soup.find("div", {'class': 'ph experiment'})
        all_attributes = line.find_all('th')
        all_att_values = line.find_all('td')
        attributes_dict = {}
        for idx, att in enumerate(all_attributes):
            try:
                attributes_dict[att.contents[0]] = all_att_values[idx].contents[0]
            except:
                continue
        currline = run
        for item in attributes:
            data = attributes_dict.get(item, 'N/A')
            currline = currline + ',' + data
        outlines.append(currline)
    with open(outfile, 'w') as outfile:
        for line in outlines:
            outfile.write('{}\n'.format(line))
    return outlines

def sort_platforms(infile, outfile1, outfile2):
    outfile1 = outfile1 + '_ion_torrent.txt'
    outfile2 = outfile2 + '_illumina.txt'
    with open(infile, 'r') as infile:
        lines = [line.strip().split(',') for line in infile]
        ion_torrent = [line[0] for line in lines if line[1] == 'Ion Torrent']
        illumina = [line[0] for line in lines if line[1] == 'Illumina']
    with open(outfile1, 'w') as out1:
        for run in ion_torrent:
            out1.write('{}\n'.format(run))
    with open(outfile2, 'w') as out2:
        for run in illumina:
            out2.write('{}\n'.format(run))
    return ion_torrent, illumina


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input SRA accession list (txt file)')
    parser.add_argument('-o', '--output', help='Prefix for all output files. Do not include extensions.')
    parser.add_argument('-a', '--attributes', nargs='+', help='List of desired attributes from each run;'
                                                              'separated by a single space')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    acc_list = sra_acc(args.input)
    pull_metadata(acc_list, args.attributes, args.output)
    sort_platforms(args.output + '_meta.csv', args.output, args.output)
