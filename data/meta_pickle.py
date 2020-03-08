import os
import sys
import argparse
import json
import pickle

'''
A simple python script to write JSON files to a dictionary in bulk and store the resulting dictionary as a pickle.
'''

def JSON_to_dict(file):
    with open(file) as infile:
        data = json.load(infile)
        return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input directory of JSON files to be pickled')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    metadata_dict = {}

    for root, dirs, files in os.walk(args.input):
        for directory in dirs:
            for root2, dirs2, files2 in os.walk(os.path.join(root, directory)):
                for name in files2:
                    if root2 == os.path.join(args.input, directory) and name[-5:] == '.json':
                        file = os.path.join(root2, name)
                        temp = JSON_to_dict(file)
                        for k, v in temp.items():
                            metadata_dict[k] = metadata_dict.get(k, []) + [v]

    metafile = open('metaPickle.p', 'ab')
    pickle.dump(metadata_dict, metafile)
    metafile.close()
