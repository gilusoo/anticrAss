import sys
import argparse

"""
This is a simple python script to generate an accession list of all of the datasets present in a particular directory. 
Temp files can be generated using the bash commands 'find ERR* > temp1.txt' and 'find SRR* > temp2.txt'
"""


def get_accs(file):
    return [x.strip().replace('_pass.fastq', '') for x in open(file).readlines()]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-1', '--temp1', help='Input list of ERR files')
    parser.add_argument('-2', '--temp2', help='Input list of SRR files')
    parser.add_argument('-o', '--output', help='Output file to write to')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)


    with open(args.output, 'w') as outfile:
        for acc in get_accs(args.temp1):
            outfile.write(f'{acc}\n')
        for acc in get_accs(args.temp2):
            outfile.write(f'{acc}\n')