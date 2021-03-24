import sys
import dendropy
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import os
import argparse
import pandas as pd

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
    parser.add_argument('--min-coverage', help='if more gap frequency than this parameter (0 to 1) then swapping score is 0', default=0,
                        type=float)
    parser.add_argument('--max-swapset', help='max number of elements in swap set', default=999, type=int)

    parser.add_argument('--cons-score-function', help='identity: percentage identity (0-1); blosum: blosum score (0-1)', default='blosum')

    parser.add_argument('--specify-cluster-center', help='must be an ID of the alignment. Only the swaps of groups from this ID will be computed', default='')
    parser.add_argument('--specify-cluster-set', help='must be a comma separated list of ID. Only the swaps of this group will be computed', default='')
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
        self.tree = self.get_tree(args)
        self.dm = self.get_distance_matrix(args)
        self.ds = self.get_sequence_dict(args)
        # seqs
        self.tree_seqs = set(self.dm.keys())
        self.seqs = set([x for x in self.ds.keys() if x in self.tree_seqs])
        # run IDs
        self.seqrun, self.seqgroup = self.get_run_IDs(args)
        # annotation
        self.annot_df = self.get_annotation(args)
        # orthologs / paralogs dictionaries
        self.orthologs, self.paralogs = self.get_matching_paralogs_dics(args)

        # output
        self.output_path = args.output
        self.output_df_counter = 0
        self.output_df_dict = {}
        self.output_df = pd.DataFrame()
        self.output_string_list = []

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

    def get_run_IDs(self, args):
        """
        handles the parameter to restrict the analysis to only one ID or a specific group of IDs
        :param args:
        :return:
        """
        seqrun = self.seqs
        if args.specify_cluster_center:
            assert args.specify_cluster_center in self.seqs
            seqrun = args.specify_cluster_center

        seqgroup = []
        if args.specify_cluster_set:
            seqgroup = args.specify_cluster_set.split(",")
            assert all([x in self.seqs for x in seqgroup])
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
                this_name_list_out = []
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

                    best_score = 0.0
                    best_report = ""
                    swap_set = set()
                    # we add two new groups for conservation measure, this_swap, other_swap
                    this_prot_swap = set()
                    other_prot_swap = set()

                    out_str = ""

                    # get the sorted list of neighbours for this ID
                    swap_list = sorted([(self.dm[name][lab], lab) for lab in list(self.dm[name].keys()) if self.dm[name][lab] <= self.MAX_DISTANCE])
                    for _, ID in swap_list:
                        if len(swap_set) >= self.MAX_SWAPSET:  # break if reach the limit of dimension of swapset
                            break
                        if self.get_protein(ID) != this_prot:  # ignore if it's not the same protein
                            continue

                        # SWAPPING LABELS
                        # remove the closest protein from the group of this protein, put in the swap group of this protein
                        this_prot_orthologs.remove(ID)
                        this_prot_swap.add(ID)

                        # remove all the others
                        for lab in self.paralogs[this_prot][ID]:
                            if lab in other_prot_orthologs:
                                other_prot_swap.add(lab)
                                other_prot_orthologs.remove(lab)
                        swap_set.add(ID)

                        # SWAP SCORE
                        this_score, report = self.get_new_score(this_prot_orthologs, other_prot_orthologs, this_prot_swap, other_prot_swap, pos, other_pos, i)
                        if this_score > best_score:
                            best_score = this_score
                            best_report = report
                            out_str = "test: {}; pos {}: {}; other_pos: {}; score: {:.2}; conservations: {:.2}-{:.2}; protein: {}; name: {}; swap_set: {}; report: {}".format(
                                this_comparison, i + 1, pos, other_pos, best_score, compare_score[0], compare_score[1], this_prot, name,
                                ",".join(swap_set), ",".join(report))

                    if best_score >= self.MIN_DIFF:  # and outer_score >= MIN_OUT_SCORE:
                        this_name_list_out.append(out_str)
                        self.output_df_dict[self.output_df_counter] = {"test": this_comparison, "pos": i + 1, "residue": pos, "pos_residue": "{}{}".format(i + 1, pos),
                                                     "other_pos": other_pos, "name": name, "score": best_score, "protein": this_prot,
                                                     "conservation": "{:.2}-{:.2}".format(compare_score[0], compare_score[1]),
                                                     "swap_set": ",".join(swap_set), "report": ",".join(best_report)}
                        self.output_df_counter += 1
                self.output_string_list.append("{}\n{}\n".format(name, "\n".join(this_name_list_out)))

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
        if pos == "-" or pos == "X":
            return 0.0, 0.0
        # find if both strings pass the min cover setting, otherwise score is zero
        non_gap_percent = lambda x: len(x.replace("-", "")) / len(x)
        if non_gap_percent(this_clone_str) < self.MIN_COVER or non_gap_percent(other_clone_str) < self.MIN_COVER:
            return 0.0, 0.0
        # calculate the conservation scores
        this_score = self.CONS_SCORE_FUNC(pos, this_clone_str)
        other_score = self.CONS_SCORE_FUNC(other_pos, other_clone_str)

        return this_score, other_score

    def get_new_score(self, this_prot, other_prot, this_prot_swap, other_prot_swap, pos, other_pos, i):
        """
        #now t7: compare the swapping pos conservation and the clone combaciante conservation, all scores are positive
        t7.5: like t7, but get the minimum of the 4 scores times 4
        :param this_clones:
        :param other_clones:
        :param ds:
        :param pos:
        :param compare_scores:
        :return:
        """

        report = ["0"] * 4
        check = pos + other_pos
        if "-" in check or "X" in check or len(check) < 2:  # if the position is a gap or one is empty return 0
            return 0.0, report

        this_prot_str = "".join([self.ds[x][i] for x in this_prot])
        other_prot_str = "".join([self.ds[x][i] for x in other_prot])
        this_prot_swap_str = "".join([self.ds[x][i] for x in this_prot_swap])
        other_prot_swap_str = "".join([self.ds[x][i] for x in other_prot_swap])

        # if the length of gapless column is zero, just return score zero
        if any([len(x.replace("-","")) == 0 for x in [this_prot_str,other_prot_str,this_prot_swap_str, other_prot_swap_str]]):
            return 0.0, report
        # find if both strings pass the min cover setting, otherwise score is zero
        non_gap_percent = lambda x: len(x.replace("-", "")) / len(x)
        if non_gap_percent(this_prot_str) < self.MIN_COVER or non_gap_percent(other_prot_str) < self.MIN_COVER:
            return 0.0, report
        if pos == other_pos:  # no swap if it's the same pos
            return 0.0, report

        # calculate the conservation scores
        this_prot_score = self.CONS_SCORE_FUNC(other_pos, this_prot_str)            # the other pos in this orthologs
        other_prot_score = self.CONS_SCORE_FUNC(pos, other_prot_str)                # this pos in the other orthologs
        this_prot_swap_score = self.CONS_SCORE_FUNC(other_pos, this_prot_swap_str)  # the other pos in this swap group
        other_prot_swap_score = self.CONS_SCORE_FUNC(pos, other_prot_swap_str)      # this pos in the other swap group

        best_score = min([this_prot_score, other_prot_score, this_prot_swap_score, other_prot_swap_score])

        report = [str(x) for x in (this_prot_score, other_prot_score, this_prot_swap_score, other_prot_swap_score)]
        return best_score, report

    def save_output(self):
        """
        save the output of the analysis
        :param args:
        :param outdf_dict:
        :return:
        """
        csv_out_file_path = self.output_path + "swapping_residues.csv"
        string_out_file_path = self.output_path + "swapping_residues.txt"
        string_out_file = open(string_out_file_path, 'w')
        string_out_file.write("".join(self.output_string_list))
        df = pd.DataFrame(self.output_df_dict)  # NEED TO ADD THE STUFF
        df.T.to_csv(csv_out_file_path, sep=",")

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
