import sys
import dendropy
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import numpy as np
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


def create_parser(argv):
    """Create a command line parser with all arguments defined"""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@')

    # General parameters
    # input
    parser.add_argument('--alg-file', help='the input alignment file')
    parser.add_argument('--tree-file', help='the input tree file')
    parser.add_argument('--annotation-file', help='the annotation file, a CSV containing as columns ID,organism,protein')

    # parameters
    parser.add_argument('--max-distance', help='the max distance to check the inversion residue cluster', default=2.8, type=float)
    parser.add_argument('--min-score-diff', help='the minimum score to report results', default=0,
                        type=float)
    parser.add_argument('--max-codon-pos', help='the last codon position to check (right limit)', default=-1,
                        type=int)
    parser.add_argument('--min-coverage', help='if more gap frequency than this parameter (0 to 1) then swapping score is 0', default=0.1,
                        type=float)
    parser.add_argument('--max-swapset', help='max number of elements in swap set', default=999, type=int)
    parser.add_argument('--pseudocounts', help='whether using the pseudocounts correction of the aa freq.', default=True)

    parser.add_argument('--cons-score-function', help='identity: percentage identity (0-1); blosum: blosum score (0-1)', default='blosum')

    parser.add_argument('--specify-cluster-center', help='must be an ID of the alignment. Only the swaps of groups from this ID will be computed', default='')
    parser.add_argument('--specify-cluster-set', help='must be a comma separated list of organisms. Only the organisms of this group will be computed', default='')
    parser.add_argument('--specify-cluster-name', help='the name of the specified set of species', default="Group X")
    parser.add_argument('--specify-root-node', help='the root node of the provided tree', default='midpoint')
    parser.add_argument('--multi-protein-approach',
                        help='The approach when more than 2 proteins (e.g. duplicates) are present. Conservative: only '+
                             'compares pairs of proteins. Inclusive: tries to match other pairs', default='Conservative')

    # output folder
    parser.add_argument('--output', help='the output folder')
    # redo
    parser.add_argument('--redo', default=False, type=bool, help='set to True if you want to remake')


    return parser


class Dirphy:

    def __init__(self, args):
        """
        init the stuff needed for analysis
        :param args:
        """
        # params
        self.MAX_DISTANCE = args.max_distance
        self.MIN_DIFF = args.min_score_diff
        self.MAX_POS = args.max_codon_pos
        self.MAX_SWAPSET = args.max_swapset

        self.MIN_COVER = args.min_coverage
        score_funcs = {
            'identity': Dirphy.find_conservation_identity,
            'blosum': Dirphy.find_conservation_blosum,
        }
        self.CONS_SCORE_FUNC = score_funcs[args.cons_score_function]
        self.MULTI_PROT_APPROACH = args.multi_protein_approach
        self.PSEUDOCOUNTS = args.pseudocounts
        # input
        self.tree = self.get_tree(args)
        self.dm = self.get_distance_matrix(args)
        self.ds = self.get_sequence_dict(args)
        # seqs
        self.tree_seqs = set(self.dm.keys())
        self.seqs = set([x for x in self.ds.keys() if x in self.tree_seqs])
        # annotation
        self.annot_df = self.get_annotation(args)
        self.org_dict = self.get_organism_dictionary()
        self.organisms = set(self.org_dict.values())
        # run IDs
        self.seqrun, self.seqgroup = self.get_run_IDs(args)
        self.groupname = args.specify_cluster_name
        # orthologs / paralogs dictionaries
        self.orthologs, self.paralogs = self.get_matching_paralogs_dics(args)

        # output
        self.output_path = args.output
        # meta output
        self.output_df_counter_meta = 0
        self.output_df_dict_meta = {}
        self.output_df_meta = pd.DataFrame()
        self.output_string_list_meta = []
        # ss out
        self.output_df_counter_ss = 0
        self.output_df_dict_ss = {}
        self.output_df_ss = pd.DataFrame()
        self.output_string_list_ss = []

        # utilities
        self.AAlist = list('ACDEFGHIKLMNPQRSTVWY')

        # pseudocounts bg and pij
        self.beta = 5
        self.BGprob = self.get_aa_background()
        self.qij_mat, self.qij_sum = self.get_qij_matrix()

        # if a group, calculate the base probability, only 2 prots
        if len([x for x in self.orthologs.keys() if x is not np.nan]) == 2:
            self.p1, self.p2 = [x for x in self.orthologs.keys() if x is not np.nan]
            self.base_prob_calc = BGprob()
            self.base_prob = self.calculate_base_prob()



    def get_sequence_dict(self, args):
        """
        get sequences from a multiple sequence alignment, check for same length requirement
        :param args:
        :return:
        """
        # open MSA
        alg = SeqIO.parse(args.alg_file, "fasta")

        # save seq in dict
        ds = {}
        for s in alg:
            ds[Dirphy.rename(s.id)] = s.seq

        # test if all have the same length
        len_first = len(list(ds.values())[0])
        assert all([len(x) == len_first for x in ds.values()])
        return ds

    def get_distance_matrix(self, args):

        tree = self.get_tree(args)  # open tree

        # distance matrix
        pdc = tree.phylogenetic_distance_matrix()
        dm = {}  # label 1 -> dictionary of label 2 -> distance

        for i, t1 in enumerate(tree.taxon_namespace):
            # if t1.label not in PRUNED_TAXA:
            this_d = {}
            for t2 in tree.taxon_namespace:
                # if t2.label not in PRUNED_TAXA:
                this_d[t2.label] = pdc(t1, t2)
            dm[t1.label] = this_d.copy()
        return dm

    def get_tree(self, args):
        try:
            tree = dendropy.Tree.get(
                path=args.tree_file,
                schema="nexus")
        except:
            tree = dendropy.Tree.get(
                path=args.tree_file,
                schema="newick")

        root_node = args.specify_root_node
        if root_node == "midpoint":
            tree.reroot_at_midpoint()
        else:
            outgroup_node = tree.find_node_with_taxon_label(Dirphy.rename(root_node))
            tree.to_outgroup_position(outgroup_node)
        return tree

    def get_annotation(self, args):
        """
        dataframe with columns: ID, organism, protein
        :param args:
        :return:
        """
        colnames = ['ID', 'organism', 'protein']
        df = pd.read_csv(args.annotation_file)
        if not all(df.columns == colnames):
            df = pd.read_csv(args.annotation_file, names=colnames)
        df['ID'] = df['ID'].apply(Dirphy.rename)
        return df

    def get_organism_dictionary(self):
        od = {}
        for i in self.annot_df.index:
            od[self.annot_df['ID'][i]] = self.annot_df['organism'][i]
        return od

    def get_run_IDs(self, args):
        """
        handles the parameter to restrict the analysis to only one ID or a specific group of IDs
        :param args:
        :return:
        """
        seqrun = self.seqs
        if args.specify_cluster_center:
            args.specify_cluster_center = Dirphy.rename(args.specify_cluster_center)
            assert args.specify_cluster_center in self.seqs, "cluster center not found in {}".format(self.seqs)
            seqrun = [args.specify_cluster_center]

        seqgroup = []
        if args.specify_cluster_set:
            seqgroup = args.specify_cluster_set.split(",")
            assert all([x in self.organisms for x in seqgroup]), "not all organisms found in the list.\norganisms: {}\nlist: {}".format(",".join(sorted(seqgroup)),",".join(sorted(self.organisms)))
        return seqrun, seqgroup

    def get_output_df(self):
        """
        if it's the first time accessing the output dataframe, create it
        """
        if not self.output_df and self.output_df_dict:
            self.output_df = pd.DataFrame(self.output_df_dict)
        return self.output_df

    # function to get the opposite position in clone combaciante and take care of multiple hits
    def get_other_pos(self, name, this_protein, other_protein, i):
        """
        returns the list of positions of all paralogs in the same organism of this protein
        :param name:
        :param protein:
        :param i:
        :return:
        """
        out = ""
        paras = list(self.paralogs[this_protein][name])
        for p in paras:
            if p in self.orthologs[other_protein]:
                out = self.ds[p][i]  # FIX HERE, for now, it will just take one of the possible (multiple) paralogs
        return out

    def get_matching_paralogs_dics(self, args):
        """
        return a dictionary of proteins orthologs, and a dictionary of proteins paralogs
        :param args:
        :param annot_df: dataframe with columns: ID, protein, organism
        :return:
        """
        # create two dictionaries of matching paralogs
        onagi_prot = {}
        matching_prot = {}
        # for every protein in the list
        prot_list = set(self.annot_df['protein'])
        for prot in prot_list:
            onagi_prot[prot] = set([Dirphy.rename(x) for x in self.annot_df.loc[self.annot_df['protein'] == prot]['ID']]) & self.seqs
            this_match_dic = {}
            for x in onagi_prot[prot]:
                this_org = str(self.annot_df.loc[self.annot_df['ID'] == x]['organism'].values[0])
                # same organism but different protein
                same_o_diff_p = (self.annot_df['protein'] != prot) & (self.annot_df['organism'] == this_org)
                this_match_dic[x] = set([Dirphy.rename(x) for x in self.annot_df.loc[same_o_diff_p]['ID']]) & self.seqs
            matching_prot[prot] = this_match_dic
        return onagi_prot, matching_prot

    # function to return clone
    def get_protein(self, name):
        protein = "no_label"
        for prot in self.orthologs:
            if name in self.orthologs[prot]:
                protein = prot
        return protein

    # function that stores amino acid bg
    def get_aa_background(self):
        """
        stores aa background in an array of len 20 with self.AAlist as labels
        now from Le S.Q., Gascuel O. Molecular Biology and Evolution. 2008 25(7):1307-20.
        :return:
        """
        # BLOSUM62 background distribution
        blosum_background_distr = [0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337,
                                   0.022355, 0.062157, 0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197,
                                   0.053287, 0.012066, 0.034155, 0.069147]
        amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        capra_bg = {amino_acids[i]:blosum_background_distr[i] for i in range(20)}
        reordered_bg = [capra_bg[AA] for AA in self.AAlist]
        return reordered_bg

    # function to get the qij matrix
    def get_qij_matrix(self):
        """
        qij is the matrix of the frequencies of co occurrence of two amino acids as calculated by the blosum formula
        qij = Pi * Pj * e^(lambda*BLOSij)
        :return:
        """

        order = [self.AAlist.index(x) for x in "A R N D C Q E G H I L K M F P S T W Y V".split()]
        O = lambda x: order[x]
        qi_raw = """
0.0215 
0.0023 0.0178 
0.0019 0.0020 0.0141 
0.0022 0.0016 0.0037 0.0213 
0.0016 0.0004 0.0004 0.0004 0.0119 
0.0019 0.0025 0.0015 0.0016 0.0003 0.0073 
0.0030 0.0027 0.0022 0.0049 0.0004 0.0035 0.0161 
0.0058 0.0017 0.0029 0.0025 0.0008 0.0014 0.0019 0.0378 
0.0011 0.0012 0.0014 0.0010 0.0002 0.0010 0.0014 0.0010 0.0093 
0.0032 0.0012 0.0010 0.0012 0.0011 0.0009 0.0012 0.0014 0.0006 0.0184 
0.0044 0.0024 0.0014 0.0015 0.0016 0.0016 0.0020 0.0021 0.0010 0.0114 0.0371 
0.0033 0.0062 0.0024 0.0024 0.0005 0.0031 0.0041 0.0025 0.0012 0.0016 0.0025 0.0161 
0.0013 0.0008 0.0005 0.0005 0.0004 0.0007 0.0007 0.0007 0.0004 0.0025 0.0049 0.0009 0.0040 
0.0016 0.0009 0.0008 0.0008 0.0005 0.0005 0.0009 0.0012 0.0008 0.0030 0.0054 0.0009 0.0012 0.0183 
0.0022 0.0010 0.0009 0.0012 0.0004 0.0008 0.0014 0.0014 0.0005 0.0010 0.0014 0.0016 0.0004 0.0005 0.0191 
0.0063 0.0023 0.0031 0.0028 0.0010 0.0019 0.0030 0.0038 0.0011 0.0017 0.0024 0.0031 0.0009 0.0012 0.0017 0.0126 
0.0037 0.0018 0.0022 0.0019 0.0009 0.0014 0.0020 0.0022 0.0007 0.0027 0.0033 0.0023 0.0010 0.0012 0.0014 0.0047 0.0125 
0.0004 0.0003 0.0002 0.0002 0.0001 0.0002 0.0003 0.0004 0.0002 0.0004 0.0007 0.0003 0.0002 0.0008 0.0001 0.0003 0.0003 0.0065 
0.0013 0.0009 0.0007 0.0006 0.0003 0.0007 0.0009 0.0008 0.0015 0.0014 0.0022 0.0010 0.0006 0.0042 0.0005 0.0010 0.0009 0.0009 0.0102 
0.0051 0.0016 0.0012 0.0013 0.0014 0.0012 0.0017 0.0018 0.0006 0.0120 0.0095 0.0019 0.0023 0.0026 0.0012 0.0024 0.0036 0.0004 0.0015 0.0196 
""".split()
        qij_mat = [[0 for i in range(20)] for j in range(20) ]
        counter = 0
        for i in range(20):
            for j in range(i+1):
                qij_mat[O(i)][O(j)] = qij_mat[O(j)][O(i)] = float(qi_raw[counter])
                counter += 1
        return qij_mat, [sum([x[i] for x in qij_mat]) for i in range(20)]

    def calculate_lengths(self):
        """
        calculates the lengths given the groupings of seqgroup
        :return:
        """

        L1,L2 = 0,0

        # L1 is the average length of the distance from root (duplication) to speciation

        tax_a1 = [x for x in self.orthologs[self.p1]]
        tax_a2 = [x for x in self.orthologs[self.p2]]
        a1 = self.tree.mrca(taxon_labels=tax_a1).distance_from_root()
        a2 = self.tree.mrca(taxon_labels=tax_a2).distance_from_root()

        L1 = np.mean([a1,a2])

        # L2 is the average length of the distance from speciation to the end
        a3 = self.tree.mrca(taxon_labels=tax_a1).distance_from_tip()
        a4 = self.tree.mrca(taxon_labels=tax_a2).distance_from_tip()

        L2 = np.mean([a3,a4])

        return L1,L2

    def calculate_base_prob(self):
        """
        calculates the base prob only if the grouping is specified
        :return:
        """

        if not self.seqgroup or len([x for x in self.orthologs if x is not np.nan]) > 2:  # for now, only specified groups of a comparison of 2
            return 0

        L1, L2 = self.calculate_lengths()
        return self.base_prob_calc.calculate_prob(L1,L2)

    def run_analysis(self):
        """
        use everything to run the analysis
        :return:
        """
        for name in self.seqrun:  # for each protein in the analysis protein list

            # skip if the protein has no label (no annotation)
            this_prot = self.get_protein(name)
            if this_prot == "no_label":
                continue

            other_proteins = [p for p in self.orthologs.keys() if p != this_prot and self.orthologs[p]]
            for other_prot in other_proteins:
                this_comparison = "{}vs{}".format(this_prot, other_prot)
                #avg_intra, avg_extra = self.calculate_average_distance(this_prot, other_prot)
                this_name_list_out_meta = []
                this_name_list_out_ss = []
                for i in range(len(self.ds[name])):  # for each position of the alignment
                    # skip position if we reached the max pos and it's not "-1"
                    if i > self.MAX_POS != -1:
                        break

                    # get the residue of this position for this protein
                    pos = self.ds[name][i]
                    # get the residue of this position for the other protein in the same organism
                    other_pos = self.get_other_pos(name, this_prot, other_prot, i)

                    # get the names of the orthologs of this protein
                    this_prot_orthologs = set(self.orthologs[this_prot])
                    # get the names of all the other proteins
                    other_prot_orthologs = set(self.orthologs[other_prot])

                    # starting conservation
                    compare_score = self.get_compare_score(this_prot_orthologs, other_prot_orthologs, pos, other_pos, i)

                    best_score_meta = 0.0
                    best_report_meta = ""
                    best_score_ss = 0.0
                    best_report_ss = ""

                    out_str_meta = ""
                    out_str_ss = ""

                    # get the sorted list of neighbours for this ID
                    swap_order = sorted([(self.dm[name][lab], lab) for lab in list(self.dm[name].keys()) if self.dm[name][lab] <= self.MAX_DISTANCE])
                    swap_list = [[swap_order[j][1] for j in range(i+1)] for i in range(min(len(swap_order),self.MAX_SWAPSET))]
                    if self.seqgroup:
                        swap_list = [[x for _,x in swap_order if self.get_protein(x) == this_prot and self.org_dict[x] in self.seqgroup]]
                    for swaps in swap_list:
                        if not all([self.get_protein(swaps[i]) == this_prot for i in range(len(swaps))]):  # ignore if it's not the same protein
                            break

                        # we add 4 new groups for conservation measure, this_ortho, other_ortho, this_swap, other_swap
                        this_prot_ortho = this_prot_orthologs.copy()
                        other_prot_ortho = other_prot_orthologs.copy()
                        this_prot_swap = set()
                        other_prot_swap = set()
                        swap_set = set()
                        # SWAPPING LABELS
                        for ID in swaps:
                            this_prot_ortho.remove(ID)
                            this_prot_swap.add(ID)

                            # remove all the others
                            for lab in self.paralogs[this_prot][ID]:
                                if lab in other_prot_ortho:
                                    other_prot_swap.add(lab)
                                    other_prot_ortho.remove(lab)
                            swap_set.add("{}".format(self.org_dict[ID]))

                        # META SCORE
                        this_score_meta, this_score_ss, report_meta, report_ss = self.get_new_score(this_prot_ortho, other_prot_ortho, this_prot_swap, other_prot_swap, i)
                        if this_score_meta > best_score_meta:
                            best_score_meta = this_score_meta
                            best_report_meta = report_meta
                            out_str_meta = "test: {}; pos {}: {}; other_pos: {}; score: {:.2}; conservations: {:.2}-{:.2}; protein: {}; name: {}; swap_set: {}; report: {}".format(
                                this_comparison, i + 1, pos, other_pos, best_score_meta, compare_score[0], compare_score[1], this_prot, name,
                                ",".join(swap_set), report_meta)
                        # SS SCORE
                        if this_score_ss > best_score_ss:
                            best_score_ss = this_score_ss
                            best_report_ss = report_ss
                            out_str_ss = "test: {}; pos {}: {}; other_pos: {}; score: {:.2}; conservations: {:.2}-{:.2}; protein: {}; name: {}; swap_set: {}; report: {}".format(
                                this_comparison, i + 1, pos, other_pos, best_score_ss, compare_score[0], compare_score[1], this_prot, name,
                                ",".join(swap_set), report_ss)
                    # update meta output
                    if best_score_meta >= self.MIN_DIFF:  # and outer_score >= MIN_OUT_SCORE:
                        this_name_list_out_meta.append(out_str_meta)
                        self.output_df_dict_meta[self.output_df_counter_meta] = {"test": this_comparison, "pos": i + 1, "residue": pos, "pos_residue": "{}{}".format(i + 1, pos),
                                                     "other_pos": other_pos, "name": name, "score": best_score_meta, "protein": this_prot,
                                                     "conservation": "{:.2}-{:.2}".format(compare_score[0], compare_score[1]),
                                                     "swap_set": ",".join(swap_set), "report": best_report_meta}
                        self.output_df_counter_meta += 1
                    # update ss output
                    if best_score_ss >= self.MIN_DIFF:  # and outer_score >= MIN_OUT_SCORE:
                        this_name_list_out_ss.append(out_str_ss)
                        self.output_df_dict_ss[self.output_df_counter_ss] = {"test": this_comparison, "pos": i + 1, "residue": pos, "pos_residue": "{}{}".format(i + 1, pos),
                                                     "other_pos": other_pos, "name": name, "score": best_score_ss, "protein": this_prot,
                                                     "conservation": "{:.2}-{:.2}".format(compare_score[0], compare_score[1]),
                                                     "swap_set": ",".join(swap_set), "report": best_report_ss}
                        self.output_df_counter_ss += 1
                self.output_string_list_meta.append("{}\n{}\n".format(name, "\n".join(this_name_list_out_meta)))
                self.output_string_list_ss.append("{}\n{}\n".format(name, "\n".join(this_name_list_out_ss)))

    def get_compare_score(self, this_prot_orthologs, other_prot_orthologs, pos, other_pos, i):
        """
        get the initial conservation score of the two clones
        :param this_clones:
        :param other_clones:
        :param ds:
        :param pos:
        :return:
        """
        this_clone_str = "".join([self.ds[x][i] for x in this_prot_orthologs])
        other_clone_str = "".join([self.ds[x][i] for x in other_prot_orthologs])
        # if the position is gap, just return score zero
        if pos == "-" or pos == "X" or other_pos == "-" or other_pos == "X":
            return 0.0, 0.0
        # find if both strings pass the min cover setting, otherwise score is zero
        non_gap_percent = lambda x: len(x.replace("-", "")) / len(x)
        if non_gap_percent(this_clone_str) < self.MIN_COVER or non_gap_percent(other_clone_str) < self.MIN_COVER:
            return 0.0, 0.0
        # calculate the conservation scores
        this_score = self.CONS_SCORE_FUNC(pos, this_clone_str)
        other_score = self.CONS_SCORE_FUNC(other_pos, other_clone_str)

        return this_score, other_score

    def get_new_score(self, this_prot, other_prot, this_prot_swap, other_prot_swap, i):
        """
        #t8: the joint probability of having a swapped residue
        :param this_clones:
        :param other_clones:
        :param ds:
        :param pos:
        :param compare_scores:
        :return:
        """
        report = ""
        this_prot_str = "".join([self.ds[x][i] for x in this_prot])
        other_prot_str = "".join([self.ds[x][i] for x in other_prot])
        this_prot_swap_str = "".join([self.ds[x][i] for x in this_prot_swap])
        other_prot_swap_str = "".join([self.ds[x][i] for x in other_prot_swap])

        # if the length of gapless column is zero, just return score zero
        if any([len(x.replace("-","")) == 0 for x in [this_prot_str,other_prot_str,this_prot_swap_str, other_prot_swap_str]]):
            return 0.0, 0.0, report, report
        # find if both strings pass the min cover setting, otherwise score is zero
        non_gap_percent = lambda x: len(x.replace("-", "")) / len(x)
        if non_gap_percent(this_prot_str) < self.MIN_COVER or non_gap_percent(other_prot_str) < self.MIN_COVER:
            return 0.0, 0.0, report, report

        # calculate the conservation scores
        best_score_meta, best_score_ss, report_meta, report_ss = self.compute_joint_probability(this_prot_str,other_prot_str,this_prot_swap_str,other_prot_swap_str)

        return best_score_meta, best_score_ss, report_meta, report_ss

    def compute_joint_probability(self, this_prot_str, other_prot_str, this_prot_swap_str, other_prot_swap_str):
        """
        function to compute the joint probabilities of having a swapped residue
        basically calculates the sum of (Ai*ci * (1-bi*di) or the same for SS

        by definition, meta is a swap between this swap (b) and other swap (c), while SS is a swap between this (a) and this swap (b)
        the result is the following matchings:
        META: a with c, b with d
        SS: a with d, b with c

           a (this_prot)
        -[
       |   b (this_swap_prot)
      -|
       |   c (other_swap_prot)
        -[
           d (other_prot)
        """
        rl = range(20)  # this should be the length of the arrays (20 aa)
        report_meta = []
        report_ss = []
        meta_score = 0
        ss_score = 0

        a = self.frequency_array(this_prot_str)
        d = self.frequency_array(other_prot_str)
        b = self.frequency_array(this_prot_swap_str)
        c = self.frequency_array(other_prot_swap_str)
        if not all([len(x) == 20 for x in [a, b, c, d]]): # full gap column...
            return 0,0,"",""

        # matching probabilities
        Pac = sum([a[i] * c[i] for i in rl])
        Pbd = sum([b[i] * d[i] for i in rl])
        Pad = sum([a[i] * d[i] for i in rl])
        Pbc = sum([b[i] * c[i] for i in rl])

        for i in range(len(a)):
            if Pac and Pbd:  # non zero matching probabilities
                Padd_meta = (a[i] * c[i] / Pac) * (1 - (b[i] * d[i] / Pbd))
                report_meta.append("{0}:{1:.2f}".format(self.AAlist[i], Padd_meta))
                meta_score += Padd_meta
            if Pad and Pbc:  # non zero matching probabilities
                Padd_ss = (a[i] * d[i] / Pad) * (1 - (b[i] * c[i] / Pbc))
                report_ss.append("{}:{:.2f}".format(self.AAlist[i], Padd_ss))
                ss_score += Padd_ss

        # with this correction, we calculate the joint probability from the conditional probability
        meta_score *= Pac * Pbd
        ss_score *= Pad * Pbc

        # report is the added prob for each AA
        report_meta = ";".join(report_meta)
        report_ss = ";".join(report_ss)
        return meta_score, ss_score, report_meta, report_ss

    def frequency_array(self, array):
        """
        return the array in single letter amino acid order of the frequency of occurrences given the list
        PSEUDOCOUNTS, when active, uses the LG matrix of aa transition and the formula of better scoring schemes..
        X AND gaps are removed!
        :param array:
        :return:
        """
        count = Counter(array)
        del(count['-'])
        del(count['X'])
        sum_arr = sum(count.values())
        if sum_arr == 0:
            return []
        freq_arr = [count.get(x, 0)/sum_arr for x in self.AAlist]
        assert round(sum(freq_arr),3) == 1, "the sum of all frequencies is not 1, {}".format(count)
        alpha = sum_arr
        beta = self.beta
        if self.PSEUDOCOUNTS:
            new_freq_arr = []
            for a in range(len(self.AAlist)):
                new_freq = freq_arr[a] * alpha / (alpha + beta)
                pseudo_sum = 0
                for i in range(len(self.AAlist)):
                    pseudo_sum += freq_arr[i] * self.qij_mat[i][a] / self.qij_sum[i]
                pseudo_sum *= beta
                new_freq += (pseudo_sum / (alpha + beta))
                new_freq_arr.append(new_freq)
            freq_arr = new_freq_arr
        return freq_arr

    def save_output(self):
        """
        save the output of the analysis
        :param args:
        :param outdf_dict:
        :return:
        """
        # meta
        csv_out_file_path = self.output_path + "inverted_residues_meta.csv"
        string_out_file_path = self.output_path + "inverted_residues_meta.txt"
        dist_plot_meta = self.output_path + "score_distribution_meta.png"
        string_out_file = open(string_out_file_path, 'w')
        string_out_file.write("".join(self.output_string_list_meta))
        df = pd.DataFrame(self.output_df_dict_meta)  # NEED TO ADD THE STUFF
        df.T.to_csv(csv_out_file_path, sep=",")
        ax = sns.displot(data=df.T, x='score')
        if self.base_prob:
            plt.axvline(self.base_prob, color='r', label="Base chance: {:02f}".format(self.base_prob))
            plt.legend()
        plt.savefig(dist_plot_meta)

        # ss
        csv_out_file_path = self.output_path + "inverted_residues_ss.csv"
        string_out_file_path = self.output_path + "inverted_residues_ss.txt"
        dist_plot_ss = self.output_path + "score_distribution_ss.png"
        string_out_file = open(string_out_file_path, 'w')
        string_out_file.write("".join(self.output_string_list_ss))
        df = pd.DataFrame(self.output_df_dict_ss)  # NEED TO ADD THE STUFF
        df.T.to_csv(csv_out_file_path, sep=",")
        ax = sns.displot(data=df.T, x='score')
        if self.base_prob:
            plt.axvline(self.base_prob, color='r', label="Base chance: {:02f}".format(self.base_prob))
            plt.legend()
        plt.savefig(dist_plot_ss)

        # logos
        if self.base_prob:
            fsize =lambda poslist: (int(np.log10(len(poslist)) * 20 + 10),
                     int(len(poslist) / 100 + 8))
            label = "{0}-Others,{0}-{2},{1}-{2},{1}-Others".format(self.p1, self.p2, self.groupname).split(",")
            seqlist = self.get_logo_sequence_list()
            # meta
            meta_logo_path = self.output_path + "meta_logo.png"
            poslist_meta = self.get_highscore_positions_meta()
            fig, axes = plt.subplots(4, 1, figsize=fsize(poslist_meta))
            for i in range(4):
                Dirphy.make_logo(axes[i], label[i], seqlist[i], poslist_meta)
            plt.tight_layout()
            plt.savefig(meta_logo_path, dpi=300)

            # ss
            ss_logo_path = self.output_path + "ss_logo.png"
            poslist_ss = self.get_highscore_positions_ss()
            fig, axes = plt.subplots(4, 1, figsize=fsize(poslist_ss))
            for i in range(4):
                Dirphy.make_logo(axes[i], label[i], seqlist[i], poslist_ss)
            plt.tight_layout()
            plt.savefig(ss_logo_path, dpi=300)

    @staticmethod
    def make_logo(ax, label, seq_list, poslist):
        import logomaker as lm
        raw_seqs = []
        for seq in seq_list:
            raw_seqs.append("".join([seq[x-1] for x in poslist]))
        counts_mat = lm.alignment_to_matrix(raw_seqs)
        normalized_mat = lm.transform_matrix(counts_mat, normalize_values=True)
        logo = lm.Logo(normalized_mat, stack_order='small_on_top',
                       font_name='Futura',
                       font_weight='ultralight',
                       color_scheme='chemistry',
                       vpad=0.1,
                       width=.9,
                       ax=ax)
        logo.style_spines(visible=False)
        logo.ax.set_yticks([0, .5, 1])
        logo.ax.set_xticks(range(len(poslist)))
        logo.ax.set_xticklabels([str(x) for x in poslist], fontsize='small')
        logo.ax.set_ylabel(label)

    def get_logo_sequence_list(self):
        """
        returns sequences for groups a, b, c, and d
        a is others-p1
        b is group1-p1
        c is group1-p2
        d is others-p2
        :return:
        """

        seqlist = []

        # a
        seqlist.append([self.ds[ID] for ID in self.orthologs[self.p1] if self.org_dict[ID] not in self.seqgroup])

        # b
        seqlist.append([self.ds[ID] for ID in self.orthologs[self.p1] if self.org_dict[ID] in self.seqgroup])

        # c
        seqlist.append([self.ds[ID] for ID in self.orthologs[self.p2] if self.org_dict[ID] in self.seqgroup])

        # d
        seqlist.append([self.ds[ID] for ID in self.orthologs[self.p2] if self.org_dict[ID] not in self.seqgroup])

        return seqlist

    def get_highscore_positions_meta(self):
        """
        get the positions for the logo, 0 based
        :return:
        """
        unsorted_poslist = []
        for i in range(self.output_df_counter_meta):
            if self.output_df_dict_meta[i]["score"] > self.base_prob:
                unsorted_poslist.append((self.output_df_dict_meta[i]["score"],self.output_df_dict_meta[i]["pos"]))
        sorted_poslist = sorted(unsorted_poslist,reverse=True)
        poslist = sorted([x for _,x in sorted_poslist[:min(50,len(sorted_poslist))]])

        return poslist

    def get_highscore_positions_ss(self):
        """
        get the positions for the logo, 0 based
        :return:
        """
        unsorted_poslist = []
        for i in range(self.output_df_counter_ss):
            if self.output_df_dict_ss[i]["score"] > self.base_prob:
                unsorted_poslist.append((self.output_df_dict_ss[i]["score"], self.output_df_dict_ss[i]["pos"]))
        sorted_poslist = sorted(unsorted_poslist, reverse=True)
        poslist = sorted([x for _, x in sorted_poslist[:min(50, len(sorted_poslist))]])

        return poslist

    @staticmethod
    def find_conservation_blosum(pos, seqlist):
        """
        i need some comments on this:
        i need a value from 0 to 1 using blosum matrix
        first i get the score by summing all substitution scores
        then i get the maximum possible score (everything conserved)
        then i get the minimum possible score (i find the lowest value in the matrix)
        then i normalize with (score - min) \ (max - min)
        """
        blosm = MatrixInfo.blosum62

        aa_seqlist = [x for x in seqlist if x != "-"]
        if len(aa_seqlist) == 1:
            return 1.0

        blos_list = []
        for x in aa_seqlist:
            try:
                blos_list.append(blosm[(pos, x)])
            except KeyError:
                blos_list.append(blosm[(x, pos)])
        blos_score = sum(blos_list)

        max_score = blosm[(pos, pos)] * len(aa_seqlist)
        min_score = min([blosm[(k)] for k in blosm.keys() if pos in k]) * len(aa_seqlist)
        # print(max_score, min_score, min([blosm[(k)] for k in blosm.keys() if pos in k ]), len(aa_seqlist)-1, (blos_score - min_score) / (max_score - min_score))
        return (blos_score - min_score) / (max_score - min_score)

    @staticmethod
    def find_conservation_identity(pos, seqlist):

        aa_seqlist = [x for x in seqlist if x != "-"]
        return aa_seqlist.count(pos) / len(aa_seqlist)

    @staticmethod    # rename so that it's the same as the phylotree names (ugly spaces ugh)
    def rename(name):
        return name.replace("_", " ").replace("|", " ")


class BGprob:
    """
    This class is used to calculate the background probability of specie-specific adaptation and residue inversion
    It's based on the simple model of a duplication, 4 branches, 2 species evolution, with equal branch lengths between
    the duplication branches (L1, 2 branches) and the speciation branches (L2, 4 branches).
    The probability of a mutation is calculated from the rate with the following formula: P = 1 - e^(-lambda).

    """
    def __init__(self):

        self.R = 1
        self.saved_probs = {}

    def get_prob(self, L1, L2):
        """
        get the probability of s-s adaptation and residue inversion from the dictionary, or if not found for the request
        L1 and L2, calculates, saves and returns it.
        :param L1:
        :param L2:
        :return:
        """

        if (L1, L2) in self.saved_probs:
            return self.saved_probs[(L1, L2)]
        else:
            return self.calculate_prob(self, L1, L2)

    def calculate_prob(self, L1, L2):
        """
        calculate the probability of s-s adaptation and residue inversion, saves it and return it
        :param L1:
        :param L2:
        :return:
        """
        prob = 0
        for i in range(64):
            pi = self.outcome_probability(i, L1=L1, L2=L2, r=self.R)
            prob = self.update_state_prob(i, pi, prob)

        self.saved_probs[(L1,L2)] = prob
        return prob

    @staticmethod
    def get_normalized_bit(val, bit_index):  # to get 0 and 1 pos
        return (val >> bit_index) & 1

    @staticmethod
    def mi(a, b, anc):  # inner match so with 1 = 1
        if a == 2:
            if b == 1:
                return 0
            else:
                return 1 / 19
        elif a == 1:
            if b == 1 and anc == 1:
                return 1
            elif b == 1 and anc == 0:
                return 1 / 19
            else:
                return 0
        elif a == 0:
            if b == 0:
                return 1
            elif b == 1:
                return 0

    @staticmethod
    def mo(a, b):  # outer match, 1 != 1
        """ this one changed, now """
        if a == 2:
            if b == 0:
                return 1 / 19
            elif b == 1:
                return (18 / 19) * (1 / 19)
            elif b == 2:
                return (18 / 19) * (18 / 19 ** 2) + (1 / 19) * (
                            1 / 19)  # if the ancestor 1 matched, then 2 matches with 1/19, otherwise
                # it matches in 18 out of 19^2 cases ( because anc is on a different state)
        elif a == 1:
            if b == 0:
                return 0
            elif b == 1:
                return 1 / 19
            else:
                return (18 / 19) * (
                            1 / 19)  # only if the ancestor of 2 didn't match with 1, they can match with prob 1/19
        else:
            if b == 2:
                return 1 / 19
            elif b == 1:
                return 0
            else:
                return 1

    def outcome_probability(self, X, L1, L2, r):
        """
        X is a string of six binary digits
        r = 1 # modifier
        L1 = rate * r # mut before species split
        L2 = rate * r # mut after species split
        """

        assert X.bit_length() <= 6, "too big number"
        L1 = L1 * r
        L2 = L2 * r
        L1 = (1 - np.e ** -L1)  # now it's a probability
        L2 = (1 - np.e ** -L2)  # now it's a probability

        p1 = [1 - L1, L1]  # each node 1 counts twice
        p2 = [1 - L2, L2]  # 0 means no mutation, with prob 1-p. 1 means mutation, with prob p

        prob = 1

        for i in range(4, 6):  # ANCESTORS
            prob *= p1[BGprob.get_normalized_bit(X, i)]
        for i in range(4):  # LEAF NODES
            prob *= p2[BGprob.get_normalized_bit(X, i)]
        # print(bin(X),prob)
        return prob

    def update_state_prob(self, X, pX, prob):
        """
        updates the outcome probability of inversion
        probability of both outer match times the min of the probabilities of inner match
        the probability of s-s and of inversion is the same and is the probability of outer match times
        """
        a0 = 0
        a1 = BGprob.get_normalized_bit(X, 5)
        a2 = BGprob.get_normalized_bit(X, 4)
        l1 = a1 + BGprob.get_normalized_bit(X, 3)
        l2 = a1 + BGprob.get_normalized_bit(X, 2)
        l3 = a2 + BGprob.get_normalized_bit(X, 1)
        l4 = a2 + BGprob.get_normalized_bit(X, 0)
        # got the leaf status

        inversion_prob = BGprob.mo(l1,l3) * BGprob.mo(l2,l4) * (1-min([BGprob.mi(l1,l2,a1),BGprob.mi(l3,l4,a2)]))  # outer match and not inner match

        prob += inversion_prob * pX
        return prob

def main(argv):

    # argparse
    try:
        parser = create_parser(argv)
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit()

    # start
    output_preparation(args)

    test = Dirphy(args)

    test.run_analysis()

    test.save_output()


def output_preparation(args):
    """
    check output path, create or redo it, save the args
    :param args:
    :return:
    """
    # make out folder
    if args.redo:
        if os.path.exists(args.output):
            for file in next(os.walk(args.output))[2]:
                os.remove(args.output + file)
            os.system('rm -r {}'.format(args.output))
    assert not os.path.exists(args.output), "output path found existing already: {}".format(args.output)
    os.makedirs(args.output)

    # save args
    with open(args.output + "args_report.txt", "w") as argfile:
        argfile.write("\n".join(["{} -> {}".format(k, v) for k, v in vars(args).items()]))


if __name__ == '__main__':
    main(sys.argv)
