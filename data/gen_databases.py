import re



def remove_tara(infile, outfile):
    with open(infile, encoding='utf8') as infile, open(outfile, 'w', encoding='utf8') as outfile:
        output = []
        for line in infile:
            line = line.strip().split('\t')
            # print(line)
            if re.findall('Tara',line[1]):
                continue
            else:
                line = '\t'.join(line)
                output.append(line)
        for line in output:
            outfile.write('{}\n'.format(line))
    return

def create_runs_dict(infile):
    ann_to_runs = {}
    with open(infile, 'r', encoding='utf8') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(',')
            ann_to_runs[line[1]] = ann_to_runs.get(line[1], []) + line[3:]
    return ann_to_runs

def pull_run_db(dict, cat_list, outfile):
    run_db = []
    for cat in cat_list:
        for run in dict.get(cat, None):
            run = run.replace('\"', '')
            run_db.append(run)
    with open(outfile, 'w') as outfile:
        for run in run_db:
            outfile.write('{}\n'.format(run))
    return

if __name__ == '__main__':

    annotations_file = 'sra_annotations.csv'
    new_annotations_file = 'SRA_annotations_notara.csv'
    condensed_file = 'SRA_notara_condensed.csv'

    remove_tara(annotations_file, new_annotations_file)

    # DB Categories
    soil = ['soil']
    gut = ['human gut']
    fresh_water = ['fresh water']
    marine = ['marine', 'marine benthic', 'marine coastal', 'marine deep', 'marine surface']

    pull_run_db(create_runs_dict(condensed_file), soil, 'soil_db.txt')
    pull_run_db(create_runs_dict(condensed_file), gut, 'gut_db.txt')
    pull_run_db(create_runs_dict(condensed_file), fresh_water, 'fresh_water_db.txt')
    pull_run_db(create_runs_dict(condensed_file), marine, 'marine_db.txt')
