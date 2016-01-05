import sys, os, re, random
from HapMixUtil import get_chrs_len
from parseVcfTruthFile import del_chrs
from parseVcfTruthFile import filter_truth_file


def get_commandline_options():
    """Get user options from command line"""

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c", "--cnv_file", dest="cnv_file", type="string", help="name of the Canvas vcf file")
    parser.add_option("-t", "--bed_truth_file", dest="bed_file",  type="string", help="name of the output bed truth file")
    parser.add_option("-o", "--output_dir", dest="output_dir",  type="string", help="name of the output directory")
    parser.add_option("-e", "--excl_file", dest="excl_file",  type="string", help="name of the bed file for filtering output")
    parser.add_option("-f", "--min_overlap_fraction", dest="ov_frac",  type="float", help="min overlap fraction for variants filtering")
    parser.add_option("-a", "--add_vars", dest="add_vars",  action="store_true", default=False, help="add variants")
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


def convert_vcf2bed(vcf_file, bed_file, hap_split=True, no_dipl=False, cn_hap_max=None):
    """Parse Canvas CNV vcf file into truth bed file"""

    with open(vcf_file, "r") as vcf:
        with open(bed_file, "w") as bed:
            bed_lines = []

            for vcf_line in vcf:

                if re.match("#", vcf_line):
                    continue

                line_el = vcf_line.strip().split("\t")
                chr_id = line_el[0]
                if not chr_id.startswith("chr"):
                    chr_id = "chr" + chr_id
                pos = line_el[2].split(":")
                start = pos[3]
                end = str(int(pos[4]) - 1)
                call = line_el[9].split(":")
                if hap_split:
                    if (random.random() < 0.5):
                        cnA, cnB = int(call[3]), int(call[2])-int(call[3])
                    else:
                        cnB, cnA = int(call[3]), int(call[2])-int(call[3])
                    if cn_hap_max is not None:
                        cnA, cnB = min(cnA, 3), min(cnB, 3)
                    if no_dipl and not (cnA==1 and cnB==1):
                        bed_lines.append("\t".join([chr_id, start, end, str(cnA), str(cnB)]))
                else:
                    cn = int(call[2])
                    if cn_hap_max is not None:
                        cn = min(cn, cn_hap_max*2)
                    bed_lines.append("\t".join([chr_id, start, end, str(cn)]))

            bed.writelines("\n".join(bed_lines))


# def add_variants(bed_file, num_var, size_var_min, size_var_max, cn_min, cn_max, hg_file):
# """Add num_var CNVs to the bed truth file (size in range size_var_min-size_var_max, cn in range cn_min-cn_max)"""
#
#         chrs_len = get_chrs_len(hg_file)
#         with open(bed_file, "r") as bed:
#             variants = bed.readlines()
#
#         with open(bed_file, "w") as bed:
#
#             to_add = num_var
#             max_it = 5000
#             added = it = 0
#
#             while added < to_add:
#
#                 it += 1
#                 if it >= max_it:
#                     print("Max iter reached for " + chr_id + "\n%d variants added, %d required" % (added, to_add))
#                     break
#
#                 chr_id = "chr" + str(random.randint(1, 22))
#                 start = random.randint(1, chrs_len[chr_id])
#                 end = min(start + random.randint(size_var_min, size_var_max), chrs_len[chr_id])
#
#                 for var in variants:
#                     (var_chr, var_start, var_end, var_cnA, var_cnB) = var.strip().split("\t")
#                     (var_start, var_end) = (int(var_start), int(var_end))
#                     if (start < var_start and end > var_start) or (start > var_start and start < var_end):
#                         print("reject " + str(start) + "-" + str(end) + ", var " + str(var_start) + "-" + str(var_end))
#                         break
#                 else:
#                     cnA = random.randint(cn_min, cn_max)
#                     cnB = random.randint(cn_min, cn_max)
#                     variants.append("\n" + "\t".join([chr_id, str(start), str(end), str(cnA), str(cnB)]))
#                     added += 1
#
#             #print("writing" + bed_file)
#             bed.writelines(variants)


def add_variants(bed_file, num_var, size_var_min, size_var_max, cn_min, cn_max, seed, hg_file):
    """Add num_var CNVs to the bed truth file (size in range size_var_min-size_var_max, cn in range cn_min-cn_max)"""

    import subprocess
    cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools complement -i %s -g %s" % (bed_file, hg_file)
    pc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    compl = pc.stdout.readlines()

    random.seed(seed)
    to_add = num_var
    max_it = 500
    added = it = 0
    variants = []
    while added < to_add:
        #print(added, to_add)
        it += 1
        if it >= max_it:
            print("Max iter reached: %d variants added, %d required" % (added, to_add))
            break
        dipl_reg = compl.pop(random.randint(1, len(compl) - 1)).strip().split("\t") # randomly select a diploid region in the bed_file
        size_var = random.randint(size_var_min, size_var_max)
        if not dipl_reg[0].startswith("chr"):
            dipl_reg[0] = "chr" + dipl_reg[0]
        if dipl_reg[0] not in ["chr" + str(chr_num) for chr_num in range(1, 23)] or (int(dipl_reg[2]) - int(dipl_reg[1])) < size_var:
            continue
        else:
            start = random.randint(int(dipl_reg[1]), int(dipl_reg[2]) - size_var)
            end = start + size_var - 1
            cnA = random.randint(cn_min, cn_max)
            cnB = random.randint(cn_min, cn_max)
            if (cnA, cnB) not in [(0, 0), (1, 1)]:
                print "adding\t" +  "\t".join([dipl_reg[0], str(start), str(end), str(cnA), str(cnB)])
                variants.append("\n" + "\t".join([dipl_reg[0], str(start), str(end), str(cnA), str(cnB)]))
                added += 1

    with open(bed_file, "a") as bed:
        bed.writelines(variants)


def main():
    """Convert Canvas output CNV file (.vcf) to truth file (.bed), optionally adding variants"""

    opts = get_commandline_options()
    hg_file = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/genomes/human.hg19.genome"

    convert_vcf2bed(opts.cnv_file, opts.bed_file, True, True, 3)
    del_chrs(opts.bed_file, ["chrX", "chrY", "chrM"], hg_file)
    if opts.add_vars:
        add_variants(opts.bed_file, 8, 1000000, 10000000, 2, 3, 99, hg_file)
    if opts.excl_file:
        filter_truth_file(opts.bed_file, opts.excl_file,opts.ov_frac)


if __name__=="__main__":
    main()


# python parseCanvasCNVFile.py -c /illumina/build/CNA/CNV/GEL_FFPE_2602_Canvas/GeL004_FF.v3/Analysis/Normal_Tumor_G1_P1.somatic.CNV.vcf -o /illumina/scratch/tmp/users/ccolombo/simulation/ -t GeL004.bed -a -e /illumina/build/CNA/CNV/Germline_Canvas/Germline/HapMixNormals/Proteus_LP6007590_LP6007591/Analysis/Normal_S1.SV.bed -f 0.3