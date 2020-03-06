import re
import argparse
import os
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Directory of BLAST results files')
    parser.add_argument('-o', '--output', help='Output csv file to write to')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    for root, dirs, files in os.walk(args.input):
        for name in files:
            if root == args.input:
                file = os.path.join(root, name)
                with open(args.output, 'a') as outfile:
                    outfile.write(f'Bin #,# of Hits,Top Hit,Score,E-value\n')
                for line in open(file):
                    if line.startswith('Query= '):
                        if curr:
                            outfile.write(f'{name},{hit_count},{top_hit},{top_score},{top_e}\n')
                        curr = line.strip().split()[1]
                        hit_count = 0
                    if line.startswith('*****'):
                        outfile.write(f'{curr},0\n')
                    if re.match(r'N[A-Z]_', line[0:2]):
                        line = line.strip().sploit()
                        if hit_count == 0:
                            top_hit = line[0]
                            top_score = line[1]
                            top_e = line[2]
                        hit_count += 1
                    if top_hit:
                        outfile.write(f'{name},{hit_count},{top_hit},{top_score},{top_e}\n')
                    else:
                        outfile.write(f'{name},0\n')