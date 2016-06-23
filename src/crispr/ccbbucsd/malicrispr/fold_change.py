"""This module exposes functions that calculate fold-changes for the Mali dual CRISPR experiments."""

# third-party libraries
import numpy

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_pseudocount_adjusted_series(data_frame, count_header):
    return data_frame[count_header].apply(lambda x: x + 1)


def add_fold_change_series(dataframe, numerator_series, numerator_name, denominator_series, denominator_name):
    fold_change = (numerator_series / denominator_series)
    fold_change_header = dataframe.make_foldchange_header(numerator_name, denominator_name)
    dataframe.add_column(fold_change, fold_change_header)
    return fold_change, fold_change_header


def add_log2fc_series(dataframe, fc_series, fc_header):
    log2_fold_change = fc_series.apply(numpy.log2)
    log2_fold_change_header = dataframe.get_log2fc_from_foldchange_header(fc_header)
    dataframe.add_column(log2_fold_change, log2_fold_change_header)


def calc_fold_changes_against_base(dataframe, denominator_header, numerators_headers):
    denominator_series = get_pseudocount_adjusted_series(dataframe, denominator_header)
    for curr_numerator_header in numerators_headers:
        # calculate fold change and log2 fold change for each gRNA pair, and add to data frame
        numerator_series = get_pseudocount_adjusted_series(dataframe, curr_numerator_header)
        fc_series, fc_header = add_fold_change_series(dataframe, numerator_series, curr_numerator_header,
                                                      denominator_series, denominator_header)
        add_log2fc_series(dataframe, fc_series, fc_header)


def add_timept_and_plasmid_fold_changes(counts_df, headers_manager, experiment_set_prefixes, earliest_timept):
    col_names = counts_df.columns.values

    # calc fold-change across timepoints for each experiment
    for curr_expt_prefix in experiment_set_prefixes:
        denominator_header = headers_manager.join_header([curr_expt_prefix, earliest_timept,
                                                          headers_manager.norm_counts_txt])
        numerators_headers = headers_manager.get_relevant_headers([curr_expt_prefix], [headers_manager.norm_counts_txt],
                                                                  col_names, [earliest_timept])
        calc_fold_changes_against_base(counts_df, denominator_header, numerators_headers)

    if headers_manager.includes_plasmid(col_names):
        # calc fold change of all days against plasmid
        denominator_header = headers_manager.join_header([headers_manager.plasmid_header,
                                                          headers_manager.norm_counts_txt])
        numerators_headers = headers_manager.get_relevant_headers(experiment_set_prefixes,
                                                                  [headers_manager.norm_counts_txt], col_names)
        calc_fold_changes_against_base(counts_df, denominator_header, numerators_headers)

