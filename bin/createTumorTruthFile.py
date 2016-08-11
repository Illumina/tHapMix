import sys, os, numpy
from HapMixUtil import *

def get_commandline_options():
    """Get user options from command line"""

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c", "--CN_distr_filename", dest="CN_distr_json", help="Distribution of copy numbers json file")
    parser.add_option("-m", "--MCC_distr_filename", dest="MCC_distr_json", help="Distribution of major chromosome counts json file")
    parser.add_option("-n", "--num_distr_filename", dest="num_distr_json", help="Distribution of number of variants json file")
    parser.add_option("-l", "--len_distr_filename", dest="len_distr_json", help="Distribution of variant lengths json file")
    parser.add_option("-o", "--output_filename", dest="output_BED", help="output tumor truth file")
    parser.add_option("-g", "--genome", dest="hg_file", type="string", help="genome file")
    parser.add_option("--seed", dest="rand_seed", help="Seed for random number generation", default=None, type=int)
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print(parser.print_help())
        sys.exit(2)
    if any([opt is None for opt in vars(options)]):
        print(parser.print_help())
        sys.exit(2)

    options.output_BED = os.path.abspath(options.output_BED)

    return options


def get_tumour_CNVs_distr(num_freq_file, CN_distr_file, MCC_cond_distr_file, len_freq_file, hg_file, seed=None):
    """Write tumor truth file"""

    import json, numpy

    with open(hg_file) as hgf:
        hg = dict([(line.strip().split("\t")[0], float(line.strip().split("\t")[1])) for line in hgf if line.strip().split("\t")[0] not in ["chrM", "chrX", "chrY"]])
    if seed:
        numpy.random.seed(seed)

    vars_dict = {}
    with open(num_freq_file) as num_file:
        num_distr = json.load(num_file)
    with open(CN_distr_file) as CN_file:
        CN_distr = json.load(CN_file)
    with open(MCC_cond_distr_file) as MCC_file:
        MCC_distr = json.load(MCC_file)
    with open(len_freq_file) as len_file:
        len_distr = json.load(len_file)

    tot = int(float(numpy.random.choice(num_distr.keys(), p=num_distr.values())))
    print tot

    num = 0
    while (num < tot):

        CN = numpy.random.choice(CN_distr.keys(), p=CN_distr.values())
        if not CN in MCC_distr:
            continue
        MCC = numpy.random.choice(MCC_distr[CN].keys(), p=MCC_distr[CN].values())
        CN_hapA, CN_hapB = int(MCC), int(CN) - int(MCC)
        if numpy.random.random() < 0.5:
            CN_hapA, CN_hapB = CN_hapB, CN_hapA

        # pick random chr
        chr_id = numpy.random.choice(hg.keys(), p=[clen/sum(hg.values()) for clen in hg.values()])
        chr_len = int(hg[chr_id])
        start = numpy.random.randint(1, chr_len)
        var_len = int(float(numpy.random.choice(len_distr.keys(), p=len_distr.values())))
        end = start + var_len

        # check if var is longer than chromosome
        if end > chr_len:
            continue

        if not chr_id in vars_dict:
            vars_dict[chr_id] = [(start, end, CN_hapA, CN_hapB)]
            num += 1
        else:
            for chr_var in vars_dict[chr_id]:
                # check if new var overlaps with any of the existing vars
                if (end >= chr_var[0] and end <= chr_var[1]) or (start >= chr_var[0] and start <= chr_var[1] or (start <= chr_var[0] and end >= chr_var[1])):
                    break
            else:
                vars_dict[chr_id].append((start, end, CN_hapA, CN_hapB))
                num += 1

    return vars_dict


def write_tumor_BED_file(output_BED, vars_dict):
    """Write tumor truth file"""

    with open(output_BED, "w") as out:
        for chr_id in vars_dict:
            for var in vars_dict[chr_id]:
                out.write("\t".join([chr_id] + [str(var_field) for var_field in var]) + "\n")


def main():

    opts = get_commandline_options()

    vars_dict = get_tumour_CNVs_distr(opts.num_distr_json, opts.CN_distr_json , opts.MCC_distr_json, opts.len_distr_json,
                                      opts.hg_file, opts.rand_seed)
    write_tumor_BED_file(opts.output_BED, vars_dict)


if __name__=="__main__":

    main()

