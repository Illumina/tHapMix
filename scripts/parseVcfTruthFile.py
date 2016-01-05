import sys, os, re, math, random


def get_commandline_options():
    """Get user options from command line"""

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c", "--cnv_file", dest="cnv_file", type="string", help="name of the input vcf file")
    parser.add_option("-t", "--bed_truth_file", dest="bed_file",  type="string", help="name of the output bed file")
    parser.add_option("-o", "--output_dir", dest="output_dir",  type="string", help="name of the output directory")
    parser.add_option("-e", "--excl_file", dest="excl_file",  type="string", help="name of the bed file for filtering output")
    parser.add_option("-f", "--min_overlap_fraction", dest="ov_frac",  type="float", help="min overlap fraction for variants filtering")
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print(parser.print_help())
        sys.exit(2)
    for opt_name, opt_value in vars(options).items():
        if opt_value is None:
            print(opt_name + " is not specified\n")
            print(parser.print_help())
            sys.exit(2)

    if not os.path.exists(options.cnv_file):
        print("\nCanvas output vcf file does not exist\n")
        sys.exit(2)
    if not os.path.exists(options.output_dir):
        print("Output directory does not exist\n")
        sys.exit(2)
    if options.excl_file is not None and options.ov_frac is None:
        print("Specify min overlap fraction for variants filtering")
        sys.exit(2)

    options.output_dir = os.path.abspath(options.output_dir)
    options.bed_file = os.path.abspath(os.path.join(options.output_dir, options.bed_file))

    return options


def convert_vcf2bed(vcf_file, bed_file):
    """Convert vcf truth file to bed truth file"""

    # input vcf file has total CN, must split into haplotype CNs in output bed file, max CN per hap is 3
    splitCN = {0: [(0, 0)],
               1: [(0, 1)],
               2: [(1, 1)], # no LOH (for now)
               3: [(1, 2)]}
    splitCN.update(split_CN(4, 20))

    with open(vcf_file, "r") as vcf:
        with open(bed_file, "w") as bed:
            bed_lines = []
            for vcf_line in vcf:
                if re.match("#", vcf_line):
                    continue
                line_el = vcf_line.strip().split("\t")
                chr = line_el[0]
                start = line_el[1]
                info = line_el[7]
                end = info.split(";")[0].split("=")[1]
                cn = info.split(";")[1].split("=")[1]
                if not (chr=="chrX" or chr=="chrY" or chr=="chrM") and not float(cn)==2:
                    cnA, cnB = random.sample(splitCN[int(math.floor(float(cn)))], 1)[0] # CN can be 0.5, 1.5, 2.5, ... in the vcf
                    if random.random() < 0.5:
                        bed_lines.append("\t".join([chr, start, end, str(cnA), str(cnB)]))
                    else:
                        bed_lines.append("\t".join([chr, start, end, str(cnB), str(cnA)]))
            bed.writelines("\n".join(bed_lines))


def split_CN(cn_min, cn_max):
    """Compute all possible combinations of haplotype CNs for each total CN"""

    splitCN = {}
    for cn_tot in range(cn_min, 7):
        splitCN[cn_tot] = [(cn_hap, cn_tot - cn_hap) for cn_hap in range(0, cn_tot) if (cn_hap < 4 and cn_tot - cn_hap < 4)]
    for cn_tot in range(7, cn_max + 1):
        splitCN[cn_tot] = splitCN[6]
    return splitCN


def del_chrs(bed_file, chr_ids, hg_file):

    with open(bed_file, "r") as bed:
        vars = bed.readlines()
        vars = [var for var in vars if var.strip().split("\t")[0] not in chr_ids]
        vars[-1] = vars[-1].rstrip()

    with open(hg_file, "r") as hg:
        for chr_line in hg:
            if len(chr_line.strip().split("\t")) != 2:
                continue
            chr_id, chr_len = chr_line.strip().split("\t")
            if chr_id in chr_ids:
                vars.append("\n" + chr_id + "\t1\t" + chr_len + "\t0\t0")

    with open(bed_file, "w") as bed:
        bed.writelines(vars)


def filter_truth_file(input_bed, excl_file, min_overlap_fraction):
    """Exclude variants from the bed file which overlap regions in the excl_file for a fraction > min_overlap_fraction"""
    import os, re

    output_bed = re.sub(".bed$", "flt.bed", input_bed)
    cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools subtract"
    cmd += " -a %s" % input_bed
    cmd += " -b %s" % excl_file
    cmd += " -N -f %f" % min_overlap_fraction
    cmd += " > %s" % output_bed
    os.system(cmd)


def main():
    """Convert vcf truth file to bed truth file"""

    opts = get_commandline_options()
    hg_file = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/genomes/human.hg19.genome"

    convert_vcf2bed(opts.cnv_file, opts.bed_file)
    del_chrs(opts.bed_file, ["chrX", "chrY", "chrM"], hg_file)
    if opts.excl_file is not None:
        filter_truth_file(opts.bed_file, opts.excl_file,opts.ov_frac)


if __name__=="__main__":
    main()


# python parseVcfTruthFile.py -c /illumina/build/CNA/CNV/CNATestData/CNATruthset/CanvasSomaticTestData/HCC2218Truth.vcf -o /illumina/scratch/tmp/users/ccolombo/simulation -t HCC2218.bed -e /illumina/build/CNA/CNV/Germline_Canvas/Germline/HapMixNormals/Proteus_LP6007590_LP6007591/Analysis/Normal_S1.SV.bed -f 0.3
# python parseVcfTruthFile.py -c /illumina/build/CNA/CNV/CNATestData/CNATruthset/CanvasSomaticTestData/HCC1187Truth.vcf -o /illumina/scratch/tmp/users/ccolombo/simulation -t HCC1187.bed -e /illumina/build/CNA/CNV/Germline_Canvas/Germline/HapMixNormals/Proteus_LP6007590_LP6007591/Analysis/Normal_S1.SV.bed -f 0.3
