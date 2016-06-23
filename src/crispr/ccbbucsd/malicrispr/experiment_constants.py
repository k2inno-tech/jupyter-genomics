"""This module exposes a class that centralizes shared constants for the Mali dual CRISPR experiments."""

# standard libraries
import copy


__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class ExperimentConstants:
    def __init__(self, dataset_name, working_dir, raw_counts_fp, plasmid_fp, count_headers, expt_set_prefixes,
                 earliest_timept, neg_control_genes, num_constructs):

        self.dataset_name = dataset_name

        self.working_dir = working_dir
        self.raw_counts_fp = raw_counts_fp
        self.plasmid_fp = plasmid_fp

        self.count_headers = copy.copy(count_headers)
        self.experiment_set_prefixes = copy.copy(expt_set_prefixes)
        self.earliest_timept = earliest_timept
        self.negative_control_genes = copy.copy(neg_control_genes)
        # per discussion with Roman 03/22/2016, I expect every raw count file will
        # include the same number of gRNA pair ids--the whole library.
        # Any gRNA pair id that was not sequenced will have a zero count.
        self.expected_num_constructs = num_constructs


# class DualCrisprDataFrame(DataFrame):
#     @property
#     def has_plasmid(self):
#         plasmid_count_headers = self.get_relevant_headers([self.plasmid_header], [self.counts_txt], None)
#         return len(plasmid_count_headers) > 0
#
#     def __init__(self, counts_dataframe, experiment_constants):
#         super().__init__(counts_dataframe)
#
#         self._constants = experiment_constants
#
#         # headers
#         self.construct_id_header = "construct_id"
#         self.plasmid_header = "plasmid"
#         self.position1_header = "gene_position1"
#         self.position2_header = "gene_position2"
#         self.gene_pair_header = "gene_pair"
#
#         # header pieces
#         self.all_expts_prefix = "all_expts"
#         self.counts_txt = "_counts_raw"
#         self.norm_counts_txt = "_counts_norm"
#         self.fold_change_txt = "_fc"
#         self.log2_fold_change_txt = "_log2fc"
#         self.log2_fitness_score_txt = "_log2fitness"
#         self.log2_pi_score_txt = "_log2pi_score"
#
#         # other miscellaneous constants
#         self.num_expected_genes = 2
#         self.pair_separator = "__"
#         self.piece_separator = "_"
#
#         self.experiment_set_prefixes = copy(experiment_constants.expt_set_prefixes)
#         self.earliest_timept = experiment_constants.earliest_timept
#         self.negative_control_genes = copy(experiment_constants.neg_control_genes)
#         # per discussion with Roman 03/22/2016, I expect every raw count file will
#         # include the same number of gRNA pair ids--the whole library.
#         # Any gRNA pair id that was not sequenced will have a zero count.
#         self.expected_num_constructs = experiment_constants.num_constructs
#
#         # header sets
#         self.general_headers = [self.construct_id_header, self.gene_pair_header,
#                                 self.position1_header, self.position2_header]
#
#         for orig_count_header in self._constants.count_headers:
#             self._rename_count_header(orig_count_header)
#
#     # def _rename_count_header(self, orig_count_header):
#     #     qualified_count_header = "{0}{1}".format(orig_count_header, self.counts_txt)
#     #     self.rename(columns={orig_count_header: qualified_count_header}, inplace=True)
#     #
#     # def get_relevant_headers(self, prefix_list, suffix_list, excluded_tokens):
#     #     if prefix_list is suffix_list is excluded_tokens is None:
#     #         raise ValueError("At least one input parameter must not be None")
#     #
#     #     headers_w_prefix = None
#     #     if prefix_list is not None:
#     #         headers_w_prefix = self.get_headers_with_prefix(prefix_list)
#     #
#     #     headers_w_prefix_and_suffix = headers_w_prefix
#     #     if suffix_list is not None:
#     #         headers_w_prefix_and_suffix = self.get_headers_with_suffix(suffix_list, headers_w_prefix)
#     #
#     #     if excluded_tokens is not None:
#     #         result = self.get_headers_without_tokens(excluded_tokens, headers_w_prefix_and_suffix)
#     #     else:
#     #         result = headers_w_prefix_and_suffix
#     #     return result
#     #
#     # def get_headers_with_prefix(self, prefix_list, header_set=None):
#     #     return self._get_headers_by_match(prefix_list, str.startswith, header_set)
#     #
#     # def get_headers_with_suffix(self, suffix_list, header_set=None):
#     #     return self._get_headers_by_match(suffix_list, str.endswith, header_set)
#     #
#     # def get_headers_without_tokens(self, token_list, header_set=None):
#     #     return self._get_headers_by_match(token_list, doesnt_contain, header_set)
#     #
#     # def join_header(self, pieces_list):
#     #     pieces_list.remove(None)
#     #     return "_".join(pieces_list)
#     #
#     # def make_foldchange_header(self,numerator_name, denominator_name):
#     #     return "{0}_vs_{1}{2}".format(numerator_name, denominator_name, self.fold_change_txt)
#     #
#     # def get_norm_from_count_header(self, header):
#     #     return header.replace(self.counts_txt, self.norm_counts_txt)
#     #
#     # def get_log2fc_from_foldchange_header(self, header):
#     #     return header.replace(self.fold_change_txt, self.log2_fold_change_txt)
#     #
#     # def get_fitness_from_log2fc_header(self, header):
#     #     return header.replace(self.log2_fold_change_txt, self.log2_fitness_score_txt)
#     #
#     # def get_pi_from_log2fc_header(self, header):
#     #     return header.replace(self.log2_fold_change_txt, self.log2_pi_score_txt)
#
#     def add_column(self, series, header):
#         add_series_to_dataframe(self, series, header)
#
#     # def get_dataframe_for_expt(self, experiment_set_prefix):
#     #     copied_df = self.copy()
#     #     all_relevant_headers = copy(self.general_headers)
#     #     expt_specific_headers = self.get_headers_with_prefix(experiment_set_prefix)
#     #     all_relevant_headers.extend(expt_specific_headers)
#     #     data_frame = copied_df[all_relevant_headers]
#     #     return data_frame
#
#     def extend_with_plasmid(self, plasmid_df):
#         self._rename_count_header(self.plasmid_header)
#         self.merge(plasmid_df, on=self.construct_id_header)
#
#     def expand_gene_info(self):
#         if not (self.position1_header in self.columns.values and self.position2_header in self.columns.values):
#             # TODO: can I blithely "apply" a method or do I have to explicitly account for "self"?
#             expanded_tuples_df = self.apply(lambda x: self._expand_info_from_construct_id(x), axis=1)
#             (self[self.position1_header], self[self.position2_header],
#              self[self.gene_pair_header]) = zip(*expanded_tuples_df)
#
#     def _expand_info_from_construct_id(self, construct_row):
#         construct_id = construct_row[self.construct_id_header]
#         construct_id_pieces = construct_id.split(self.pair_separator)
#
#         if len(construct_id_pieces) != self.num_expected_genes:
#             raise ValueError("Construct {0} split into {1} pieces rather than {2}".format(construct_id,
#                                                                                           len(construct_id_pieces),
#                                                                                           self.num_expected_genes))
#
#         results = []
#         for index in range(0, self.num_expected_genes):
#             gene_pieces = construct_id_pieces[index].split(
#                 self.piece_separator)
#             gene = gene_pieces[0]
#
#             # TODO: Double-check with customer if this is correct:
#             # if the "gene" name is for one of the non-targeting
#             # controls, compress to the control's name
#             for control_gene in self.negative_control_genes:
#                 if gene.startswith(control_gene):
#                     gene = control_gene
#
#             results.append(gene)
#
#         # Ensure that gene pairs are treated as combinations not
#         # permutations (order not important) by sorting gene list
#         # alphabetically
#         gene_pair = self.pair_separator.join(sorted(results))
#
#         results.append(gene_pair)  # UNsorted
#         return tuple(results)
#
#     # def _get_headers_by_match(self, items_to_match, match_func, header_set=None):
#     #     header_set = header_set if not header_set is None else self.columns.values
#     #     result = []
#     #
#     #     for curr_item in items_to_match:
#     #         matching_headers = [x for x in header_set if match_func(x, curr_item)]
#     #         result.extend(matching_headers)
#     #
#     #     return result
