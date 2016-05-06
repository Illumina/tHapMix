import os
import sys
import subprocess
import re
from optparse import OptionParser
import numpy
from parseCanvasCNVFile import convert_vcf2bed
from HapMixUtil import add_BED_complement
from HapMixUtil import get_clone_ploidy
from getVarPerc import create_var_perc_file
from getVarPerc import create_het_cn_file


def init_param_truth_file(truth_file, out_dir, hg_file):

    # Calculate tumor purity and ploidy
    # overall_ploidy = get_clone_ploidy(sorted_compl_cnv_file, hg_file, ["chrX", "chrY", "chrM"])

    # Complement truth file with diploid regions
    new_truth_file = os.path.join(out_dir, "cntot_" + os.path.basename(truth_file))
    with open(truth_file, "r") as tf:
        tf_lines = tf.readlines()
        if len(tf_lines[1].strip().split("\t")) == 5:
            new_tf_lines = []
            for tf_line in tf_lines:
                (chr_id, start, end, cnA, cnB) = tf_line.strip().split("\t")
                new_tf_lines.append("\t".join([chr_id, start, end, str(int(cnA)+int(cnB))]))
        with open(new_truth_file, "w") as new_tf:
            new_tf.writelines("\n".join(new_tf_lines))
        truth_file = new_truth_file
    cnv_file = add_BED_complement(truth_file, hg_file, sort=True, out_dir=out_dir, hap_split=False)
    #cnv_file = re.sub(".bed$", "_cnv.bed", out_file)
    subprocess.Popen("rm %s" % re.sub("_sorted.bed", "_compl.bed", cnv_file), shell=True)
    if os.path.exists(new_truth_file):
        subprocess.Popen("rm %s" % new_truth_file, shell=True)

    return cnv_file


def init_param_canvas_file(canvas_cnv_file, out_file):

    # Read tumor purity and tumor ploidy
    tumor_purity = 0
    overall_ploidy = 0
    with open(canvas_cnv_file, "r") as vcf:
        for line in vcf:
            # if line.startswith("##OverallPloidy"):
            #     overall_ploidy = float(line.strip().split("=")[1])
            if line.startswith("##EstimatedTumorPurity"):
                tumor_purity = float(line.strip().split("=")[1])
                break
            # if (tumor_purity!=0 and overall_ploidy!=0):
            #     break

    # Read estimanted haplotype coverage
    hap_coverage = 0
    if os.path.exists(os.path.join(os.path.dirname(canvas_cnv_file), "/Logging/SomaticCNV-G1_P1.stdout.txt")):
        with open(os.path.join(os.path.dirname(canvas_cnv_file), "/Logging/SomaticCNV-G1_P1.stdout.txt"), "r") as log_file:
            for log_line in log_file:
                if log_line.startswith(">>> Refined model"):
                    hap_coverage = float(log_line.strip().split(",").split(" ")[2])/2
                    print "Haplotype coverage: %f" % hap_coverage
                    break

    # Convert file format for CNV calls from vcf to bed
    cnv_bed_file = re.sub(".bed$", "_cnv.bed", out_file)
    convert_vcf2bed(canvas_cnv_file, cnv_bed_file, False)

    return cnv_bed_file, tumor_purity, hap_coverage


def create_ev_cov_files(cnv_file, cnv_type, part_file, out_dir, out_filename, tumor_purity, excl_file, save_tmp, add_perc, add_shift, hg_file):

    out_file = os.path.join(out_dir, out_filename)
    out_file_bin = re.sub(".bed$", "_bin.bed", out_file)
    hap_coverage = 0
    if cnv_type == "canvas":
        cnv_file, tumor_purity, hap_coverage = init_param_canvas_file(cnv_file, out_file)
    elif cnv_type == "truth":
        cnv_file = init_param_truth_file(cnv_file, out_dir, hg_file)
    rm_files = [cnv_file]
    print "Tumor purity: %f" % tumor_purity

    # Exclude regions
    if excl_file is not None:
            if not os.path.exists(excl_file):
                print "File for excluded regions %s does not exist\n" % excl_file
            else:
                cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools subtract"
                cmd += " -a %s" % cnv_file
                cmd += " -b %s" % excl_file
                cnv_file = run_bedtools_cmd(cmd, cnv_file, "_excl")
                rm_files.append(cnv_file)

    # Extract partitioned file
    extr_part_file = re.sub(".bed$", "_partitioned.bed", out_file)
    pc = subprocess.Popen("gunzip -c %s > %s" % (part_file, extr_part_file), shell=True)
    pc.wait()

    # Intersect variants file with bin file
    cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools intersect"
    cmd += " -a %s" % cnv_file
    cmd += " -b %s" % extr_part_file
    cmd += " -wa -wb"
    intersect_file = run_bedtools_cmd(cmd, out_file, "_int")
    rm_files.extend([intersect_file, extr_part_file])

    with open(intersect_file, "r") as ifile:
        int_lines = ifile.readlines()
        if hap_coverage==0:
            hap_coverage = numpy.median([float(int_line.strip().split("\t")[7]) for int_line in int_lines if float(int_line.strip().split("\t")[3])==2])/2
            print "Haplotype coverage: %f (median)" % hap_coverage
        wlines = ["#chr\tstart\tend\tCN\tobservedCoverage\texpectedCoverage"]
        for int_line in int_lines:
            (chr_segm, start_segm, end_segm, cn_segm, chr_bin, start_bin, end_bin, cov_bin, segm_bin) = int_line.strip().split("\t")
            expected_coverage = round((float(cn_segm)*tumor_purity + 2*(1-tumor_purity)) * hap_coverage, 3)
            wlines.append("\t".join([chr_bin, start_bin, end_bin, cn_segm, cov_bin, str(expected_coverage)]))
    with open(out_file_bin, "w") as ofile:
        ofile.writelines("\n".join(wlines))

    # # Sort chromosome intersect file
    # cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools sort"
    # cmd += " -i %s" % intersect_file
    # sorted_file = run_bedtools_cmd(cmd, out_file, "_sort")
    # rm_files.append(sorted_file)
    # Sort chromosome intersect file
    cmd = "sort -k 1,1 -k2,2n %s" % intersect_file
    sorted_file = run_bedtools_cmd(cmd, out_file, "_sort")
    rm_files.append(sorted_file)

    # Merge intersect file (median observed coverage)
    cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools merge"
    cmd += " -i %s" % intersect_file
    cmd += " -d -1"
    cmd += " -c 4,8"
    cmd += " -o mean,median"
    merged_file = run_bedtools_cmd(cmd, intersect_file, "_merg")
    rm_files.append(merged_file)

    # Calculate expected coverage for each segment
    with open(merged_file, "r") as mfile:
        wlines = ["#chr\tstart\tend\tCN\tobservedCoverage\texpectedCoverage"]
        for cnv_line in mfile.readlines():
            (chrid, start, end, cn, cov) = cnv_line.strip().split("\t")
            expected_coverage = round((float(cn)*tumor_purity + 2*(1-tumor_purity)) * hap_coverage, 3)
            wlines.append(cnv_line.strip() + "\t" + str(expected_coverage))
    with open(out_file, "w") as ofile:
        ofile.writelines("\n".join(wlines))

    # Calculate expected coverage for each segment with shifted
    cn_shift_suffix = []
    if add_shift:
        for cn_shift in [-2, -1, 1, 2]:
            with open(merged_file, "r") as mfile:
                wlines = ["#chr\tstart\tend\tCN\tobservedCoverage\texpectedCovearge"]
                for cnv_line in mfile.readlines():
                    (chrid, start, end, cn, cov) = cnv_line.strip().split("\t")
                    expected_coverage = round(((float(cn) + cn_shift)*tumor_purity + 2*(1-tumor_purity)) * hap_coverage, 3)
                    wlines.append(cnv_line.strip() + "\t" + str(expected_coverage))
                cn_shift_suffix.append("_CNshift_" + str(cn_shift))
                with open(re.sub(".bed$", cn_shift_suffix[-1] + ".bed", out_file), "w") as ofile:
                    ofile.writelines("\n".join(wlines))

    # Remove temporary files
    if not save_tmp:
        for rmf in rm_files:
            print "rm %s" % os.path.join(out_dir, rmf)
            subprocess.Popen("rm %s" % os.path.join(out_dir, rmf), shell=True)

    print "Output bed file written"

    if add_perc:
        if os.path.basename(out_file)[:-4].endswith("_truth"):
            sim_dir = "/illumina/scratch/tmp/users/ccolombo/simulation/" + os.path.basename(out_file)[4:-10]
        else:
            sim_dir = "/illumina/scratch/tmp/users/ccolombo/simulation/" + os.path.basename(out_file)[4:-4]
        if os.path.exists(sim_dir):

            # For each variant: clonal percentage, coverage of overlapping segments
            perc_file = create_var_perc_file(sim_dir, None, True)
            cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools intersect"
            cmd += " -a %s" % perc_file
            cmd += " -b %s" % out_file
            cmd += " -wa -wb"
            run_bedtools_cmd(cmd, out_file, "_perc")

            # For the whole genome: clonal percentage, coverage of overlapping bins
            cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools complement"
            cmd += " -i %s" % perc_file
            cmd += " -g %s" % hg_file
            pc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            compl = pc.stdout.readlines()
            with open(perc_file, "r") as pbed:
                lines = pbed.readlines()
            for line in compl:
                lines.append("\n" + line.rstrip() + "\t" + str(1) + "\t" + str(1) + "\t" + str(1.0))
            with open(perc_file, "w") as pbed:
                pbed.writelines(lines)
            cmd = "/illumina/thirdparty/bedtools/bedtools2-2.22.1/bin/bedtools intersect"
            cmd += " -a %s" % perc_file
            cmd += " -b %s" % out_file_bin
            cmd += " -wa -wb"
            run_bedtools_cmd(cmd, out_file_bin, "_perc_all")
            print("rm %s" % perc_file)
            subprocess.Popen("rm %s" % perc_file, shell=True)
            print "Output bed file with variants percentages written"

        else:
            print "Cannot find simulation directory"



def run_bedtools_cmd(cmd, base_out_file, suffix):

    print cmd
    out_file = re.sub(".bed$", suffix + ".bed", base_out_file)
    with open(out_file, "w") as wbed:
        pc = subprocess.Popen(cmd, stdout=wbed, shell=True)
        pc.wait()
    return out_file


def main():

    parser = OptionParser()
    parser.add_option("-i", "--CNV_file", dest="cnv_file", type="string", help="copy number variations file")
    parser.add_option("-c", "--CNV_type", dest="cnv_type", type="choice", choices=["canvas", "truth"], help="Evaluating coverage for Canvas output or for truth file? ['canvas'|'truth']")
    parser.add_option("-d", "--partitioned_file", dest="part_file", type="string", help="partitioned file (coverage for each bin)")
    parser.add_option("-f", "--out_file", dest="out_file", type="string", help="output base filename")
    parser.add_option("-o", "--out_dir", dest="out_dir", type="string", help="output directory")
    parser.add_option("-t", "--tum_pur", dest="tum_purity", type="float", help="tumor purity (specify only for truth file)")
    parser.add_option("-e", "--excl_file", dest="excl_file", type="string", help="optional bed file for excluded regions")
    parser.add_option("-p", "--add_perc", dest="add_perc",  action="store_true", default=False, help="output complete file with variants percentages")
    parser.add_option("-s", "--add_shift", dest="add_shift",  action="store_true", default=False, help="output file for shifted CN Canvas call")
    (options, args) = parser.parse_args()

    if not os.path.exists(options.cnv_file):
        if os.path.exists(options.cnv_file + ".gz"):
            print options.cnv_file + ".gz"
            pc = subprocess.Popen("gunzip %s" % options.cnv_file + ".gz", shell=True)
            pc.wait()
    if options.cnv_file is None or not os.path.exists(options.cnv_file):
        parser.error("No file or invalid file for copy number variations: %s" % options.cnv_file)
    if not (options.cnv_file.endswith(".vcf") or options.cnv_file.endswith(".bed")):
        sys.exit("invalid file format for cnv file\n")
    if options.part_file is None or not os.path.exists(options.part_file):
        parser.error("No file or invalid file for partition (coverage for each bin) : %s" % options.part_file)
    if options.out_file is None:
        options.out_file = "cov_" + options.part_file.split("/")[-4] + ".bed"
        if options.cnv_file.endswith(".bed"):
            options.out_file = "cov_" + options.part_file.split("/")[-4] + "_truth.bed"
        print "Using default base output name: " + options.out_file
    if options.out_dir is None or not os.path.exists(options.out_dir):
        options.out_dir = "/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov/cov_" + options.part_file.split("/")[-4]
        if not os.path.exists(options.out_dir):
            os.mkdir(options.out_dir)
        print "Using default output directory: " + options.out_dir
    if options.cnv_type=="canvas" and options.tum_purity:
        print "Tumor purity for canvas calls will be read from Canvas output, supplied value will be disregarded"
    if options.cnv_type=="truth" and not options.tum_purity:
        parser.error("Tumor purity must be specified")

    hg_file = "/home/ccolombo/filtered_human.hg19.genome"
    save_tmp = False

    create_ev_cov_files(options.cnv_file, options.cnv_type, options.part_file, options.out_dir, options.out_file, options.tum_purity, options.excl_file, save_tmp, options.add_perc, options.add_shift, hg_file)


if __name__ == "__main__":
    main()


# python /home/ccolombo/HapMix/getCanvasCovBed.py -c /illumina/scratch/tmp/users/ccolombo/Canvas/simGeL004_n1_m-1_v-1/Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf -d /illumina/scratch/tmp/users/ccolombo/Canvas/simGeL004_n1_m-1_v-1/Analysis/TempCNV_G1_P1/G1_P1.partitioned
# python /home/ccolombo/HapMix/getCanvasCovBed.py -c /illumina/scratch/tmp/users/ccolombo/simulation/GeL004.bed -d /illumina/scratch/tmp/users/ccolombo/Canvas/simGeL004_n1_m-1_v-1/Analysis/TempCNV_G1_P1/G1_P1.partitioned


