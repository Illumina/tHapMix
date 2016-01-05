import sys, os, numpy
from HapMixUtil import *

def get_commandline_options():
    """Get user options from command line"""

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-n", "--num_variants_min", dest="num_var_min", type="int", help="minimum number of variants in each chromosome")
    parser.add_option("-m", "--num_variants_max", dest="num_var_max", type="int", help="maximum number of variants in each chromosome")
    parser.add_option("-s", "--size_variants_min", dest="size_var_min", type="int", help="minimum size of variants")
    parser.add_option("-t", "--size_variants_max", dest="size_var_max", type="int", help="maximum size of variants")
    parser.add_option("-c", "--cn_min", dest="cn_min", type="int", help="minimum copy number")
    parser.add_option("-d", "--cn_max", dest="cn_max", type="int", help="maximum copy number")
    parser.add_option("-o", "--output_filename", dest="output_BED", help="name of the output tumor truth file")
    parser.add_option("-g", "--genome", dest="hg_file", type="string", help="genome file")
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print(parser.print_help())
        sys.exit(2)
    if any([opt is None for opt in vars(options)]):
        print(parser.print_help())
        sys.exit(2)

    options.output_BED = os.path.abspath(options.output_BED)

    return options


def write_tumor_BED_file(output_BED, num_var_min, num_var_max, size_var_min, size_var_max, cn_min, cn_max, hg_file):
    """Write tumor truth file"""

    with open(hg_file) as hgf:
        hg = dict([line.strip().split("\t") for line in hgf])
    print hg

    print output_BED
    with open(output_BED, "w") as out:
        lines = []
        for chr_id in ["chr" + str(num) for num in range(1,23)]:
            len_chr = int(hg[chr_id])
            print len_chr
            num_var = numpy.random.randint(num_var_min, num_var_max+1)
            start = sorted(numpy.random.randint(1, len_chr+1, num_var))
            next = [x-1 for x in start[1:]] + [len_chr]
            for s, n in zip(start, next):
                size = numpy.random.randint(size_var_min, size_var_max+1)
                cnA = numpy.random.randint(cn_min, cn_max+1)
                if cnA == 0:
                    cnB = numpy.random.randint(max(1, cn_min), cn_max+1)
                else:
                    cnB = numpy.random.randint(cn_min, cn_max+1)
                lines.append("\t".join([chr_id, str(s), str(min(s+size, n)), str(cnA), str(cnB)]))# + "\n")
        #lines[-1] = lines[-1].rstrip()
        out.writelines("\n".join(lines))


def main():

    opts = get_commandline_options()
    write_tumor_BED_file(opts.output_BED, opts.num_var_min, opts.num_var_max, opts.size_var_min, opts.size_var_max,
                         opts.cn_min, opts.cn_max, opts.hg_file)


if __name__=="__main__":

    main()

