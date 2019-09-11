import re


class Blast_Reader:
    """
    This is a class containing functions to help parse through results files created by NCBI's BLAST for the command line.
    """
    def __init__(self, filename):
        self.name = filename
        self.lines = [x.strip() for x in open(filename)]

    def hit_summary(self):
        summary = [line.strip().split() for line in open(self.name) if re.match(r"[A-Z]{3}[0-9]{4}", line[0:8])]
        return summary

    def hit_scores(self, summary):
        prot_to_e = {}
        for line in summary:
            prot_to_e[line[0]] = line[-1]
        return prot_to_e

    def prot_hits(self, summary):
        prot_to_e = {}
        for line in summary:
            prot_to_e[line[0]] = ' '.join(line[1:-2])
        return prot_to_e

    def write_hits(self, e_dict, bin, outfile):
        with open(outfile, 'w') as outfile:
            outfile.write('Bin,Protein,e-Value\n')
            for k, v in e_dict:
                if len(re.findall(r"e-0[3456789]", v) or re.findall(r"e-[123456789]{2}", v)) >= 1:
                    outfile.write(f'{bin},{k},{v}\n')
        return

# br = Blast_Reader('test5.out')
# # print(br.hit_summary())
# # for k, v in (br.hit_scores(br.hit_summary()).items()):
# #     print(k, v)
# br.write_hits(br.hit_scores(br.hit_summary()).items(), 'run3_bin5', 'testy.test')
