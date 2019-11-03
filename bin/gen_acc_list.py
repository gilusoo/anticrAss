"""
This script allows for the generation of a random accession list to be used within the anticrAss pipeline. SRA runs
from a large database file are randomly selected. A web scrape is then performed to ensure correct sequencing platform.
Runs that fit the necessary criteria are then written to an output file for future use.
"""

import sys
import random
import requests
from bs4 import BeautifulSoup
import argparse

def db_to_list(infile):
    with open(infile, 'r') as infile:
        acc_list = [line.strip() for line in infile]
    return acc_list

def rand_run(acc_list):
    return acc_list[random.randint(0, len(acc_list) - 1)]

def scrape(run):
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
    return attributes_dict

def check_and_append(run, att_dict, platform, layout='Single'):
    if att_dict['Platform'].startswith(platform):
        if att_dict['Layout'] == layout.upper():
            return run
        else:
            return None
    else:
        return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input SRA database file (txt file)')
    parser.add_argument('-o', '--output', help='Output filename for accession list')
    parser.add_argument('-n', '--num', type=int, help='Number of runs to include in accession list')
    parser.add_argument('-p', '--platform', help='Desired seq platform (Illumina or Ion)')
    parser.add_argument('-l', '--layout', help='Desired library layout (Single or Paired)')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    acc_list = db_to_list(args.input)
    count = 0
    with open(args.output, 'a') as outfile:
        while count < args.num:
            run = rand_run(acc_list)
            if check_and_append(run, scrape(run), args.platform):
                outfile.write('{}\n'.format(run))
                count += 1
            else:
                continue

