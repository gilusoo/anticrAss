import sys
import os
import argparse

'''
A simple python script to split FASTA files so that each resulting file contains only the sequencing data from a
single header. This will be used to split genome database files into individual files to be analyzed using FastANI
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input FASTA file to be split')
    parser.add_argument('-o', '--output', help='Output directory to write files to')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    path = args.output

    with open(args.input) as infile:
        curr = ''
        seq = []
        count = 0
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if curr != '':
                    with open(os.path.join(path, f'{str(count)}.fasta'), 'w') as outfile:
                        outfile.write(f'{curr}\n')
                        for item in seq:
                            outfile.write(f'{item}\n')
                    count += 1
                    curr = line
                seq = []
                curr = line
            else:
                seq.append(line)
        with open(os.path.join(path, f'{str(count)}.fasta'), 'w') as outfile:
            outfile.write(f'{curr}\n')
            for item in seq:
                outfile.write(f'{item}\n')
