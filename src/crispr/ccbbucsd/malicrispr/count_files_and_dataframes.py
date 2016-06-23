# standard libraries
import os

# third-party libraries
import pandas

# ccbb libraries
from ccbbucsd.utilities.analysis_run_prefixes import strip_run_prefix

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_counts_df(counts_fp, run_prefix):
    curr_counts_df = pandas.read_table(counts_fp, comment="#")  # TODO: remove hardcode
    orig_count_header = curr_counts_df.columns.values[-1]  # last column
    _, revised_count_header = os.path.split(orig_count_header)
    revised_count_header = strip_run_prefix(revised_count_header, run_prefix)
    curr_counts_df.rename(columns={orig_count_header: revised_count_header}, inplace=True)
    return revised_count_header, curr_counts_df