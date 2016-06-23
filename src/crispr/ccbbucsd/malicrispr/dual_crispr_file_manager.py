# standard libraries
import copy

# third-party libraries
import pandas

# ccbb libraries
from ccbbucsd.utilities.files_and_paths import make_file_path

# project-specific libraries
from ccbbucsd.malicrispr.construct_splitter import expand_gene_info

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def rename_count_header(dataframe, orig_count_header, count_text):
    if not orig_count_header.endswith(count_text):
        qualified_count_header = "{0}{1}".format(orig_count_header, count_text)
        dataframe.rename(columns={orig_count_header: qualified_count_header}, inplace=True)


def rename_count_headers(dataframe, count_headers, count_text):
    for orig_count_header in count_headers:
        rename_count_header(dataframe, orig_count_header, count_text)


def get_dataframe_for_header_prefix(input_df, header_prefix, header_mgr):
    copied_df = input_df.copy()
    all_relevant_headers = copy.copy(header_mgr.general_headers)
    expt_specific_headers = header_mgr.get_headers_with_prefix(header_prefix)
    all_relevant_headers.extend(expt_specific_headers)
    data_frame = copied_df[all_relevant_headers]
    return data_frame


class DualCrisprFileManager:
    def __init__(self, headers_mgr, experiment_constants):
        self.headers_mgr = headers_mgr
        self.constants = experiment_constants

        # file name pieces
        self.norm_suffix = "_norm.csv"
        self.log2_fold_change_suffix = "_log2_fc.csv"
        self.pi_suffix = "_pi_scores.csv"
        self.by_gene_pair_suffix = "_by_gene_pair.csv"

    def get_counts_dataframe(self, raw_counts_fp, plasmid_fp):
        # TODO: make sure I pass in correct header set
        raw_counts_df = pandas.read_table(raw_counts_fp, names=self.constants.count_headers, skiprows=1)
        rename_count_headers(raw_counts_df, self.constants.count_headers, self.headers_mgr.count_text)
        expand_gene_info(raw_counts_df, self.headers_mgr, raw_counts_df.columns.values,
                         self.constants.negative_control_genes)

        if not self.constants.plasmid_fp is None:
            plasmid_df = pandas.read_table(plasmid_fp,
                                           names=[self.headers_mgr.construct_id_header, self.headers_mgr.plasmid_header])
            rename_count_header(plasmid_df, self.headers_mgr.plasmid_header, self.headers_mgr.count_text)
            raw_counts_df.merge(plasmid_df, on=self.headers_mgr.construct_id_header)

        return raw_counts_df

    def write_normalized_to_file(self, df_to_save, subset_name=None):
        return self._write_df_to_file_by_suffix(self.norm_suffix, df_to_save, subset_name)

    def write_foldchanges_to_file(self, df_to_save, subset_name=None):
        return self._write_df_to_file_by_suffix(self.log2_fold_change_suffix, df_to_save, subset_name)

    def _get_file_path_by_suffix(self, suffix, subset_name):
        subset_name = self.constants.dataset_name if subset_name is None else subset_name
        return make_file_path(self.constants.working_dir, subset_name, suffix)

    def _write_df_to_file_by_suffix(self, suffix, df_to_save, subset_name):
        output_path = self._get_file_path_by_suffix(suffix, subset_name)
        df_to_save.to_csv(output_path, index=False)
        return output_path
