import sys
import argparse
import random
from bs4 import BeautifulSoup
import requests
from time import sleep


def rand_run(acc_list):
    '''
    Produces a random accession from a provided list or database file.
    :param acc_list: Input file with one accession per line
    :return: A single random accession (str)
    '''
    acc_list = [x.strip() for x in acc_list.readlines()]
    return acc_list[random.randint(0, len(acc_list) - 1)]


def scrape(run):
    '''
    Perform a web scrape for checking metadata.
    :param run:
    :return:
    '''
    req = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={}'.format(run))
    soup = BeautifulSoup(req.text, 'lxml')
    line = soup.find("div", {'class': 'ph experiment'})
    try:
        all_attributes = line.find_all('th')
        all_att_values = line.find_all('td')
    except:
        return False
    attributes_dict = {}
    for idx, att in enumerate(all_attributes):
        try:
            attributes_dict[att.contents[0]] = all_att_values[idx].contents[0]
        except:
            continue
    return attributes_dict


def check_and_append(acc, att_dict, platform, layout='Single'):
    '''
    Checks that metadata matches criteria and adds appropriate runs to output accession list
    :param acc:
    :param att_dict:
    :param platform:
    :param layout:
    :return:
    '''
    if att_dict['Platform'].startswith(platform) and att_dict['Layout'] == layout.upper():
        return acc
    else:
        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input SRA accession list or database file with one accession per line')
    parser.add_argument('-o', '--output', help='Output txt file to write accession list to')
    parser.add_argument('-n', '--num', type=int, help='Number of runs to include in accession list')
    parser.add_argument('-p', '--platform', help='Desired seq platform (Illumina or Ion)')
    parser.add_argument('-l', '--layout', help='Desired library layout (Single or Paired)')
    parser.add_argument('-r', '--random', help='Use to have the acc list generated in random order from input file')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    acc_list = [x.strip() for x in open(args.input).readlines()]
    count = 0
    used = []
    with open(args.output, 'a') as outfile:
        if args.random:
            if count < args.num:
                run = rand_run(acc_list)
                if run not in used and check_and_append(run, scrape(run), args.platform):
                    outfile.write('{}\n'.format(run))
                    count += 1
                    sleep(1)
        else:
            for run in acc_list:
                if count < args.num and scrape(run) and check_and_append(run, scrape(run), args.platform, args.layout):
                    outfile.write('{}\n'.format(run))
                    count += 1
                    sleep(1)
                else:
                    continue
