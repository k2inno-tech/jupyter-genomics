"""This module exposes functions that perform DESeq-style median normalization for the Mali dual CRISPR experiments."""

# third-party libraries
import numpy
import scipy.stats.mstats

# ccbb libraries
from ccbbucsd.utilities.pandas_utils import add_series_to_dataframe

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def add_normalization_across_expts(counts_df, headers_manager, expected_num_constructs):
    geo_mean_header = "geo_mean_across_expts"
    size_factor_txt = "_size_factor"

    # calculate geometric mean across all counts cols for each
    # gRNA pair id and add as new column
    count_columns = headers_manager.get_raw_count_headers(counts_df.columns.values)
    raw_counts_across_counts_cols = counts_df[count_columns]
    geo_means = raw_counts_across_counts_cols.apply(scipy.stats.mstats.gmean, axis=1)  # NB: 1 = apply to each row
    add_series_to_dataframe(counts_df, geo_means, geo_mean_header)

    # for each day
    for curr_count_header in count_columns:
        # calculate raw count value / geometric mean for all
        # gRNA pair ids in a count col, then remove infinite and NaN values
        raw_div_by_geo_mean = (counts_df[curr_count_header] / counts_df[geo_mean_header])
        temp = raw_div_by_geo_mean.replace([numpy.inf, -numpy.inf], numpy.nan)
        raw_div_by_geo_mean = temp.dropna()

        # calculate size factor for count col as median of above
        size_factor = raw_div_by_geo_mean.median(axis=0)  # NB: 0 = apply to each column

        # add new column holding size factor for current count col (same for all gRNA pair ids)
        size_factor_vector = [size_factor for x in range(0, expected_num_constructs)]
        size_factor_header = "{0}{1}".format(curr_count_header, size_factor_txt)
        add_series_to_dataframe(counts_df, size_factor_vector, size_factor_header)

        # vector-calculate normalized count value for all gRNA pair ids
        # in this count column as raw count value for gRNA pair id in this count column
        # divided by size factor for count column, add as new column
        norm_counts = (counts_df[curr_count_header] / counts_df[size_factor_header])
        norm_header = headers_manager.get_norm_from_count_header(curr_count_header)
        add_series_to_dataframe(counts_df, norm_counts, norm_header)





