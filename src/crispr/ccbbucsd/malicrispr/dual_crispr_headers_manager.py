__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def doesnt_contain(main_str, substr):
    return substr in main_str


class DualCrisprHeadersManager:
    def __init__(self):
        # headers
        self.construct_id_header = "construct_id"
        self.plasmid_header = "plasmid"
        self.position1_header = "gene_position1"
        self.position2_header = "gene_position2"
        self.gene_pair_header = "gene_pair"

        # header pieces
        self.all_expts_prefix = "all_expts"
        self.counts_txt = "_counts_raw"
        self.norm_counts_txt = "_counts_norm"
        self.fold_change_txt = "_fc"
        self.log2_fold_change_txt = "_log2fc"
        self.log2_fitness_score_txt = "_log2fitness"
        self.log2_pi_score_txt = "_log2pi_score"

        # header sets
        self.general_headers = [self.construct_id_header, self.gene_pair_header,
                                self.position1_header, self.position2_header]

    def includes_plasmid(self, columns_list):
        plasmid_count_headers = self.get_relevant_headers([self.plasmid_header], [self.counts_txt], columns_list,
                                                          None)
        return len(plasmid_count_headers) > 0

    def get_raw_count_headers(self, header_set):
        return self.get_headers_with_suffix([self.counts_txt], header_set)

    def get_relevant_headers(self, prefix_list, suffix_list, columns_list, excluded_tokens):
        if prefix_list is suffix_list is excluded_tokens is None:
            raise ValueError("At least one input parameter must not be None")

        headers_w_prefix = None
        if prefix_list is not None:
            headers_w_prefix = self.get_headers_with_prefix(prefix_list, columns_list)

        headers_w_prefix_and_suffix = headers_w_prefix
        if suffix_list is not None:
            headers_w_prefix_and_suffix = self.get_headers_with_suffix(suffix_list, headers_w_prefix)

        if excluded_tokens is not None:
            result = self.get_headers_without_tokens(excluded_tokens, headers_w_prefix_and_suffix)
        else:
            result = headers_w_prefix_and_suffix
        return result

    def get_headers_with_prefix(self, prefix_list, header_set):
        return self._get_headers_by_match(prefix_list, str.startswith, header_set)

    def get_headers_with_suffix(self, suffix_list, header_set):
        return self._get_headers_by_match(suffix_list, str.endswith, header_set)

    def get_headers_without_tokens(self, token_list, header_set):
        return self._get_headers_by_match(token_list, doesnt_contain, header_set)

    def join_header(self, pieces_list):
        pieces_list.remove(None)
        return "_".join(pieces_list)

    def make_foldchange_header(self, numerator_name, denominator_name):
        return "{0}_vs_{1}{2}".format(numerator_name, denominator_name, self.fold_change_txt)

    def get_norm_from_count_header(self, header):
        return header.replace(self.counts_txt, self.norm_counts_txt)

    def get_log2fc_from_foldchange_header(self, header):
        return header.replace(self.fold_change_txt, self.log2_fold_change_txt)

    def get_fitness_from_log2fc_header(self, header):
        return header.replace(self.log2_fold_change_txt, self.log2_fitness_score_txt)

    def get_pi_from_log2fc_header(self, header):
        return header.replace(self.log2_fold_change_txt, self.log2_pi_score_txt)

    def _get_headers_by_match(self, items_to_match, match_func, header_set):
        result = []
        for curr_item in items_to_match:
            matching_headers = [x for x in header_set if match_func(x, curr_item)]
            result.extend(matching_headers)

        return result
