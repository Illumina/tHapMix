#!/illumina/thirdparty/python/Python-2.7.3/bin/python
import sys
import os
import json
import re


def calc_tot_var_perc(sim_dir):
    """For a single simulation: for each variant in the truth file, calculate total percentage of sample with it (variant heterogeneity)"""

    vars = {}
    with open("%s/sim_param.json" % sim_dir) as sim_param:
        clones_perc = json.load(sim_param)["clones_perc"]

    for clone, perc in clones_perc.items():
        clone_bed_fn = [file for file in os.listdir(sim_dir) if file.endswith("_clone{0}.bed".format(str(clone)))][0]
        with open(os.path.join(sim_dir, clone_bed_fn), "r") as clone_bed:
            for line in clone_bed:
                var = line.strip()
                if var not in vars.keys():
                    vars[var] = perc
                else:
                    vars[var] += perc
    return vars


def create_var_perc_file(sim_dir, out_file, no_dipl):
    """For a single simulation: write bed file with variants percentages"""

    vars =  calc_tot_var_perc(sim_dir)
    if out_file is None:
        out_file = os.path.join(sim_dir, "var_perc_" + os.path.basename(sim_dir) + ".bed")
    with open(out_file, "w") as perc_bed:
        for var, perc in vars.items():
            (chr_id, start, end, cnA, cnB) = var.split("\t")
            if (no_dipl and not (cnA=="1" and cnB=="1")) or not no_dipl:
                perc_bed.write(var + "\t" + str(perc) + "\n")
    return out_file


def create_het_cn_file(sim_dir, out_file, no_dipl, hap_split=False, round=False):
    """Calculate variant heterogeneous CN from variant percentage file (hetCN_var = perc_var*CN_var + (1-perc_var)*2"""

    import math
    vars =  calc_tot_var_perc(sim_dir)
    if out_file is None:
        out_file = os.path.join(sim_dir, "het_cn_" + os.path.basename(sim_dir) + ".bed")
    with open(out_file, "w") as cn_bed:
        for var, perc in vars.items():
            (chr_id, start, end, cnA, cnB) = var.split("\t")
            if (no_dipl and not (cnA=="1" and cnB=="1")) or not no_dipl:
                cnA, cnB = float(cnA)*float(perc) + 2*(1 - float(perc)), float(cnB)*float(perc) + 2*(1 - float(perc))
                if round:
                    cnA, cnB = int(math.ceil(cnA)), int(math.ceil(cnB))
                if hap_split:
                    cn = str(cnA) + "\t" + str(cnB)
                else:
                    cn = str(cnA + cnB)
                cn_bed.write("\t".join([chr_id, start, end, cn]) + "\n")
    return out_file



def main(argv):

    sim_id = argv[0]
    main_dir = "/illumina/scratch/tmp/users/ccolombo/simulation/"
    sim_dirs = [main_dir + dir for dir in os.listdir(main_dir) if os.path.basename(dir).startswith(sim_id)]

    for sim_dir in sim_dirs:
        print sim_dir
        print create_var_perc_file(sim_dir, None, True)
        print create_het_cn_file(sim_dir, None, 0.8, True)


if __name__ == "__main__":
    main(sys.argv[1:])