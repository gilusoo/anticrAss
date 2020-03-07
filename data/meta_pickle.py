import os
import sys
import argparse
import json
import pickle

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
        for name in files:
            if root == args.input and name[-5:] == '.json':
                file = os.path.join(root, name)
                temp = JSON_to_dict(file)
                metadata_dict = {**metadata_dict, **temp}

    metafile = open('metaPickle.p', 'ab')
    pickle.dump(metadata_dict, metafile)
    metafile.close()
