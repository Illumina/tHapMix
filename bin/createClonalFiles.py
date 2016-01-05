import sys, os, random, re, math, json


def get_commandline_options():
    """Get user options from command line"""

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--tree_fixed_random", dest="tree_format", type="choice", choices=["fixed", "random_binary", "random_single_level"], help="tree format: random or fixed")
    parser.add_option("-t", "--tum_truth_file", dest="tum_truth_file", type="string", help="truth file")
    parser.add_option("-c", "--clones", dest="clones", type="string", action="append", help="names of the clones")
    parser.add_option("-u", "--tree_structure_file", dest="tree_file", type="string", help="tree structure file for fixed tree format")
    parser.add_option("-e", "--random_tree_seed", dest="tree_seed", type="int", help="seed to initialize generation of random tree (optional)")
    parser.add_option("-v", "--perc_heterogeneous_variants", dest="var_het", type="float", help="percentage of heterogeneous variants")
    parser.add_option("-m", "--mutate_somatic_variants", dest="mut_som_snv", action="store_true", default=False, help="add clonal somatic variants")
    parser.add_option("-s", "--somatic_snv_file", dest="som_snv_file", type="string", default="", help="name of the somatic snvs file")
    parser.add_option("-o", "--output_dir", dest="output_dir", type="string", help="name of the output directory")
    parser.add_option("-g", "--genome", dest="hg_file", type="string", help="genome file")
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print(parser.print_help())
        sys.exit(2)
    if any([opt is None for opt in vars(options)]):
        print(parser.print_help())
        sys.exit(2)

    if not os.path.exists(options.tum_truth_file):
        print "Tumor truth file does not exist\n"
        sys.exit(2)
    if options.tree_file and not os.path.exists(options.tree_file):
        print "Tree structure file does not exist\n"
        sys.exit(2)
    if options.tree_format=="fixed" and not options.tree_file:
        print "Tree structure file for fixed tree format not supplied\n"
        sys.exit(2)
    if options.tree_format=="random_binary" and options.tree_file:
        print "Tree structure file specified for random tree: file will be disregarded\n"
    if options.tree_format=="fixed" and options.tree_seed:
        print "Seed specified for fixed tree: value will be disregarded\n"
    if options.mut_som_snv and not os.path.exists(options.som_snv_file):
        print "Somatic snv file does not exist\n"
        sys.exit(2)

    options.output_BED = os.path.join(options.output_dir, os.path.basename(options.tum_truth_file)[:-4] + "_clone")
    if options.mut_som_snv:
        options.output_VCF = os.path.join(options.output_dir, os.path.basename(options.som_snv_file)[:-4] + "_clone")

    return options


def build_random_binary_clonal_tree(tree, seed=None):
    """Build a random binary clonal evolutionary tree from a list of clones"""

    import random
    if seed: random.seed(seed)
    while len(tree) > 1:
        ix = random.sample(range(0, len(tree)), 2)
        new = (tree[ix[0]], tree[ix[1]])
        del tree[max(ix)]
        del tree[min(ix)]
        tree.append(new)
    return tree[0]


def pre_order(binary_tree_node, prev_node_file, tot_vars, num_ss, base_output_file_name, format):
    """Visit binary clonal evolutionary tree and subsample variants"""

    if isinstance(binary_tree_node, (list, tuple)) and len(binary_tree_node) > 1:
        print binary_tree_node
        pre_order(binary_tree_node[0], subsample_rows(binary_tree_node[0], prev_node_file, tot_vars, num_ss, base_output_file_name, format), tot_vars, num_ss, base_output_file_name, format)
        pre_order(binary_tree_node[1], subsample_rows(binary_tree_node[1], prev_node_file, tot_vars, num_ss, base_output_file_name, format), tot_vars, num_ss, base_output_file_name, format)


def visit_non_binary_tree(tree, root_file, tot_vars, num_ss, base_output_file_name, format):
    """Visit non-binary, 1-level clonal evolutionary tree and subsample variants"""

    for tree_node in tree:
        subsample_rows(tree_node, root_file, tot_vars, num_ss, base_output_file_name, format)


def subsample_rows(node, prev_file_name, tot_rows, num_ss, base_output_file_name, format):
    """Subsample variants from tumor truth file or from somatic snv file at each node of the clonal evolutionary tree"""

    if prev_file_name is "": # root node
        curr_rows = []
    else:
        with open(prev_file_name, "r") as prev_file:
            curr_rows = [r.strip() for r in prev_file.readlines()]

    print len(tot_rows)
    print num_ss
    ind_vars = random.sample(range(len(tot_rows)), num_ss)
    for i in sorted(ind_vars, reverse=True):
        curr_rows.append(tot_rows.pop(i))

    output_file_name = base_output_file_name + "".join(re.findall("\w", str(node))) + format
    with open(output_file_name, "w") as curr_file:
        curr_file.writelines("\n".join(curr_rows))

    return output_file_name



def write_normal_BED_file(hg_file, base_output_BED_name):
    """Write BED file for normal cells"""

    with open(hg_file) as genome:
        with open(re.sub("_clone", "_norm", base_output_BED_name) + ".bed", "w") as bed:
            bed_lines = []
            for gen_line in genome:
                (chr, end) = gen_line.rstrip().split("\t")
                bed_lines.append("\t".join([chr, "0", end, "1", "1"]))
            bed.writelines("\n".join(bed_lines))


def write_clonal_files(clones, tree, var_het, tum_truth_file, output_BED, mut_som_snv, som_snv_file, output_VCF, output_dir, hg_file):
    """Write truth files (and optionally snv files) for each clone"""

    num_clones = len(clones)

    # Monoclonal scenario
    if num_clones == 1:
        with open(tum_truth_file, "r") as input_bed:
            blines = input_bed.readlines()
        with open(output_BED + clones[0] + ".bed", "w") as output_bed:
            output_bed.writelines(blines)
        if mut_som_snv:
            with open(som_snv_file, "r") as input_vcf:
                vlines = input_vcf.readlines()
            with open(output_VCF + clones[0] + ".vcf", "w") as output_vcf:
                output_vcf.writelines(vlines)
        write_normal_BED_file(hg_file, output_BED)
        return

    # Tree structures:
    if len(tree) == 2: # multi-level binary tree
        tree_structure = "binary"
        num_edges = 2*(num_clones - 1)
    elif len(tree) > 2: # single-level non-binary tree
        tree_structure = "single_level"
        num_edges = num_clones

    # Calculate percentage of variants for each branch
    if var_het >= 0:
        perc_ss_root = 1 - var_het/100 # perc common variants
        perc_ss = var_het/100/num_edges # perc variants to subsample at each node
    else: # split evenly variants also at root
        perc_ss_root = perc_ss = 1/(float(num_edges) + 1)
    print(perc_ss_root)
    print(perc_ss)

    # Bed files for cnvs (vcf files for snvs): input tumor and output clonal
    files_in = [tum_truth_file]
    files_out = [output_BED]
    if mut_som_snv:
        files_in.append(som_snv_file)
        files_out.append(output_VCF)

    # Create clonal files
    for (file_in, file_out) in zip(files_in, files_out):

        format = os.path.splitext(file_in)[1]
        print("Creating clonal " + format + " files...")

        with open(file_in, "r") as input:
            if format == ".vcf":
                tot_rows = [r.strip() for r in input.readlines() if not r.strip()[0] == "#"]
            else:
                tot_rows = [r.strip() for r in input.readlines()]
        num_vars_tot = len(tot_rows)

        root_file = subsample_rows(tree, "", tot_rows, int(math.floor(perc_ss_root*num_vars_tot)), file_out, format)
        if num_clones > 1:
            if tree_structure == "binary":
                pre_order(tree, root_file, tot_rows, int(math.floor(perc_ss*num_vars_tot)), file_out, format)
            elif tree_structure == "single_level":
                visit_non_binary_tree(tree, root_file, tot_rows, int(math.floor(perc_ss*num_vars_tot)), file_out, format)

        # Remove old and temporary files
        pattern = os.path.basename(file_out) + "(.*)" + format + "$"
        all_clone_files = [f for f in os.listdir(output_dir) if re.search(pattern, f) is not None]
        rm_files = [f for f in all_clone_files if re.search(pattern, f).group(1) not in clones]
        os.system("rm " + " ".join([os.path.abspath(os.path.join(output_dir, f)) for f in rm_files]))

    # Add normal
    write_normal_BED_file(hg_file, output_BED)



def main():

    opt = get_commandline_options()

    # Build tree for different tree formats
    if opt.tree_format == "fixed":
        tree = json.load(open(opt.tree_file))
    elif opt.tree_format == "random_binary":
        tree =  build_random_binary_clonal_tree(opt.clones[:], opt.tree_seed)
    elif opt.tree_format == "random_single_level":
        tree =  opt.clones

    # Write clonal truth files
    write_clonal_files(opt.clones, tree, opt.var_het, opt.tum_truth_file, opt.output_BED,
                       opt.mut_som_snv, opt.som_snv_file, opt.output_VCF, opt.output_dir, opt.hg_file)



if __name__=="__main__":

    main()
