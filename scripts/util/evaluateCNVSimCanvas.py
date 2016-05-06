import os, sys, re
from optparse import OptionParser
from getVarPerc import create_var_perc_file
from getVarPerc import create_het_cn_file
from HapMixUtil import add_BED_complement


def evaluate_CNV(res_dir, truth_file, excl_file, out_file, roi_file=None):
    """Run EvaluateCNV on Canvas results in res_dir and print results in files in the EvaluateCNV directory"""

    if not os.path.exists("%s/Analysis/Normal_Tumor_G1_P1.report.pdf" % res_dir): # Canvas is not finished
        print("\nNo Canvas results in %s\n\n" % res_dir)
    else:
        if not os.path.exists("%s/Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf" % res_dir):
            os.system("gunzip %s/Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf.gz" % res_dir)
        ev_cmd = "/illumina/development/tools/EvaluateCNV/1.0.1.0/EvaluateCNV %s" % truth_file
        ev_cmd += " %s/Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf" % res_dir
        #ev_cmd += " /illumina/scratch/tmp/users/sivakhno/GermlineSimulation/excluded.bed"
        ev_cmd += " %s" % excl_file
        ev_cmd += " %s" % out_file
        print ev_cmd
        os.system(ev_cmd)


def get_sim_params(sim_dir, truth_file, perc_file, out_file, hg_file):

        # Write file with simulation parameters
        # FIXME configuration file
        sim_pars = sim_dir.split("/")[-1].split("_")
        truth = sim_pars[0][3:]
        sim_pars_str = "\n".join(["%s=%s" % (sim_par[0], sim_par[1:]) for sim_par in sim_pars[1:4]])
        if len(sim_pars)==5:
            tree="sl"
        else:
            tree="bin"
        num_var = size_var_tot = size_var_het = size_gen_tot = 0
        with open(perc_file, "r") as pf:
            for p_line in pf:
                num_var += 1
                (chrid, start, end, cnA, cnB, perc) = p_line.strip().split("\t")
                size_var_tot += int(end) - int(start)
                if perc!="1.0":
                    size_var_het += int(end) - int(start)
        with open(hg_file, "r") as gf:
            for g_line in gf:
                (chr_id, chr_len) = g_line.strip().split("\t")
                size_gen_tot += int(chr_len)
        with open(out_file, "w") as ef:
            ef.writelines("\n".join(["truth=%s" % truth,
                                     sim_pars_str,
                                     "tree=%s" % tree,
                                     "numVariants=%d" % num_var,
                                     "sizeVariantsTot=%d" % size_var_tot,
                                     "sizeVariantsHet=%d" % size_var_het,
                                     "percHet=%.2f\n" % (size_var_het/float(size_gen_tot))]))


def main():
    """Evaluate Canvas results on simulated data"""

    parser = OptionParser()
    parser.add_option("-s", "--canvas_sim_id", dest="canvas_sim_id", type="string", help="simulation id (find Canvas results directories strating with sim_id")
    parser.add_option("-d", "--canvas_sim_dir", dest="canvas_sim_dir", type="string", help="Canvas results directory")
    (options, args) = parser.parse_args()

    if options.canvas_sim_dir is None and options.canvas_sim_id is None:
        parser.error("Specify Canvas results simulation id or directory")
    if options.canvas_sim_dir is not None and options.canvas_sim_id is not None:
        parser.error("Canvas results simulation id and directory are mutually exclusive options")
    if options.canvas_sim_dir is not None:
        res_dirs = [options.canvas_sim_dir]
    if options.canvas_sim_id is not None:
        canvas_dir = "/illumina/scratch/tmp/users/ccolombo/Canvas/"
        res_dirs = [canvas_dir + cdir for cdir in os.listdir(canvas_dir) if os.path.isdir((canvas_dir + cdir)) and cdir.startswith(options.canvas_sim_id)]
        if len(res_dirs) == 0:
            print "No Canvas results directories found"
    hg_file = "/home/ccolombo/filtered_human.hg19.genome"

    for res_dir in res_dirs:

        print res_dir

        sim_dir = re.sub("/Canvas/", "/simulation/", res_dir)
        res_id = re.sub("sim", "ev", res_dir.split("/")[-1])

        if not os.path.exists(res_dir):
            print("\nCanvas results directory %s does not exist\n\n" % res_dir)
            continue
        if not os.path.exists(res_dir):
            print("\nSimulation directory %s does not exist\n\n" % sim_dir)
            continue
        if os.path.basename(res_dir).startswith("simNorm") or res_dir.endswith("_purity100"):
            continue

        out_dir = os.path.join("/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCNV", res_id)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        truth_file = create_het_cn_file(sim_dir, os.path.join(out_dir, "het_truth_file.bed"), no_dipl=True, hap_split=True, round=True)
        #truth_file = "/illumina/scratch/tmp/users/ccolombo/simulation/HCC2218flt.bed"
        truth_file = add_BED_complement(truth_file, hg_file, sort=False, out_dir=out_dir, hap_split=True)
        perc_file = create_var_perc_file(sim_dir, os.path.join(out_dir, "var_perc.bed"), no_dipl=True)

        # Write simulation parameters to file
        out_file = os.path.join(out_dir, res_id + "_par.txt")
        if not os.path.exists(out_file):
            get_sim_params(sim_dir, truth_file, perc_file, out_file, hg_file)

        # Run EvaluateCNV
        out_file = os.path.join(out_dir, res_id + ".txt")
        excl_file = "/illumina/scratch/tmp/users/ccolombo/evaluation/sim_filter.bed"
        if not os.path.exists(out_file):
        #if True:
            evaluate_CNV(res_dir, truth_file, excl_file, out_file)

        # Run EvaluateCNV only on heterogeneous variants
        out_file = os.path.join(out_dir, res_id + "_onlyhet.txt")
        #if True:
        if not os.path.exists(out_file):
            with open(excl_file, "r") as ef:
                excl_vars = ef.readlines()
            with open(perc_file, "r") as pf:
                for line in pf:
                    (chr_id, start, end, cnA, cnB, perc) = line.strip().split("\t")
                    if float(perc) >= 0.8 and chr_id not in ["chrX", "chrY", "chrM"]:
                        excl_vars.append("\t".join([chr_id, start, end]) + "\n")
            os.system("rm %s" % (perc_file))
            excl_file = os.path.join(out_dir, "filter_onlyhet.bed")
            with open(excl_file, "w") as wf:
                wf.writelines(excl_vars)

            evaluate_CNV(res_dir, truth_file, excl_file, out_file)


if __name__ == "__main__":
    main()
