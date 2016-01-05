"""
Utility functions for the HapMix workflow
"""


def validate_tree_structure(tree_format, num_clones):
    """Verify if tree structure is admissable"""

    import sys

   # fixed structure
    if isinstance(tree_format, (list, tuple)):
        # single-level
        if len(tree_format) > 2:
            if not all([isinstance(tree_el, (str, unicode)) for tree_el in tree_format]):
                sys.exit("Invalid tree format!\nAllowed formats: binary tree, single-level tree\n")
            clones = tree_format
        # binary
        elif len(tree_format) == 2:
            num_children = []
            set_num_children(tree_format, num_children)
            if not all([num == 2 for num in num_children]):
                sys.exit("Invalid tree format!\nAllowed formats: binary tree, single-level tree\n")
            clones = get_clone_names_binary(tree_format)
            # check correctness
        if len(clones) != num_clones:
            sys.exit("Number of evolutionary tree leaves does not match number of clones\n")
        if len(clones) != len(set(clones)):
            sys.exit("Duplicated clone name\n")

    elif tree_format not in ["random_binary", "random_single_level"]:
        sys.exit("Invalid tree format\n")


def get_clones_names(tree_format, num_clones):
    """Return list of clones in avolutionary tree"""

    import sys
    clones = []

    # random tree
    if tree_format in ["random_binary", "random_single_level"]:
        labels = list("ABCDEFGHIJKLMNOPQRSTWXYZ") # assign default clone names
        if num_clones==2:
            clones = ["MajCl", "A"]
        else:
            clones = ["MajCl"] + labels[:(num_clones-1)]
    # fixed structure
    elif isinstance(tree_format, (list, tuple)):
        # monoclonal
        if len(tree_format) == 1:
            clones = tree_format
        else:
            # single-level
            if len(tree_format) > 2:
                clones = tree_format
            # binary
            elif len(tree_format) == 2:
               clones = get_clone_names_binary(tree_format)

    return clones


def get_clone_names_binary(binary_tree):
    """Return list of clones present in a binary evolutionary tree (leaves)"""

    tree_list = binary_tree[:]
    for i, x in enumerate(tree_list):
        while isinstance(tree_list[i], list):
            tree_list[i:i+1] = tree_list[i]
    return tree_list


def set_num_children(binary_tree_node, num_children):
    """Visit binary clonal evolutionary tree and return list of number of children for each non-leaf node """

    if isinstance(binary_tree_node, (list, tuple)) and len(binary_tree_node) > 1:
        num_children.append(len(binary_tree_node))
        set_num_children(binary_tree_node[0], num_children)
        set_num_children(binary_tree_node[1], num_children)


def add_BED_complement(truth_BED_file, hg_file, sort=True, out_dir=None, hap_split=True):
    """Add complement diploid regions to truth bed files for CNVs"""

    import subprocess, re, os
    if out_dir is None:
        out_dir = os.path.dirname(truth_BED_file)

    # Add complement regions to truth_BED_file, set to diploid
    cmd = "bedtools complement -i %s -g %s" % (truth_BED_file, hg_file)
    pc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    compl = pc.stdout.readlines()
    with open(truth_BED_file, "r") as rbed:
        lines = [line.strip() for line in rbed.readlines()]
        for line in compl:
            if hap_split:
                lines.append(line.rstrip() + "\t" + str(1) + "\t" + str(1))
            else:
                lines.append(line.rstrip() + "\t" + str(2))
    out_bed = os.path.join(out_dir, re.sub(".bed$","_compl.bed", os.path.basename(truth_BED_file)))
    with open(out_bed, "w") as wbed:
        wbed.writelines("\n".join(lines))

    # Sort output bed file
    if sort:
        cmd = "bedtools sort -i %s" % out_bed
        pc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        sortedbed = pc.stdout.readlines()
        out_bed = os.path.join(out_dir, re.sub(".bed$","_sorted.bed", os.path.basename(truth_BED_file)))
        with open(out_bed, "w") as wbed:
            wbed.writelines(sortedbed)

    return out_bed



def get_clone_ploidy(clone_truth_file, hg_file, chr_rm=[]):
    """Compute overall ploidy of sample from truth bed file"""

    clone_ploidy = 0

    with open(clone_truth_file, "r") as tfile:
        for var_line in tfile:
            (chr_id, start, end, cnA, cnB) = var_line.strip().split("\t")
            if chr_id not in chr_rm:
                clone_ploidy += (int(end) - int(start)) * (float(cnA) + float(cnB))

    gen = 0
    with open(hg_file, "r") as gfile:
        for chrl in gfile:
            (chr_id, len_chr) = chrl.strip().split("\t")
            if chr_id not in chr_rm:
                gen += int(len_chr)

    clone_ploidy = clone_ploidy/gen

    return clone_ploidy



def get_chrs_len(hg_file):
    """Read genome file and return chromosome lengths"""

    chrs = {}
    gen_file = open(hg_file, "r")
    for line in gen_file:
        if len(line.strip()) > 0:
            (chr_id, chr_len) = line.strip().split("\t")
            chrs[chr_id] = int(chr_len)
    return chrs



