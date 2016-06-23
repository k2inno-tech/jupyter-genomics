# third-party libraries
import matplotlib.pyplot
import numpy
import pandas

# ccbb libraries
from ccbbucsd.utilities.analysis_run_prefixes import strip_run_prefix
from ccbbucsd.utilities.files_and_paths import build_multipart_fp, get_file_name_pieces, get_filepaths_by_prefix_and_suffix

# project-specific libraries
from ccbbucsd.malicrispr.count_files_and_dataframes import get_counts_df

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"

DEFAULT_PSEUDOCOUNT = 1


def get_boxplot_suffix():
    return "boxplots.png"


def make_log2_series(input_series, pseudocount_val):
    revised_series = input_series + pseudocount_val
    log2_series = revised_series.apply(numpy.log2)
    nan_log2_series = log2_series.replace([numpy.inf, -numpy.inf], numpy.nan)
    return nan_log2_series.dropna().reset_index(drop=True)
    # note that .reset_index(drop=True) is necessary as matplotlib boxplot function (perhaps among others)
    # throws an error if the input series doesn't include an item with index 0--which can be the case if
    # that first item was NaN and was dropped, and series wasn't reindexed.
    

def show_and_save_histogram(output_fp, title, count_data):
    matplotlib.pyplot.figure(figsize=(20,20))
    matplotlib.pyplot.hist(count_data)
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.xlabel("log2(raw counts)")
    matplotlib.pyplot.ylabel("Frequency")
    matplotlib.pyplot.savefig(output_fp)
    matplotlib.pyplot.show()


def show_and_save_boxplot(output_fp, title, samples_names, samples_data, rotation_val=0):
    fig = matplotlib.pyplot.figure(1, figsize=(20,20))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(samples_data)
    ax.set_xticklabels(samples_names, rotation=rotation_val)   
    ax.set_xlabel("samples")
    ax.set_ylabel("log2(raw counts)")

    matplotlib.pyplot.title(title)
    fig.savefig(output_fp, bbox_inches='tight')
    matplotlib.pyplot.show()


def plot_raw_counts(input_dir, input_run_prefix, counts_suffix, output_dir, output_run_prefix, boxplot_suffix):
    counts_fps_for_run = get_filepaths_by_prefix_and_suffix(input_dir, input_run_prefix, counts_suffix)
    
    for curr_counts_fp in counts_fps_for_run:
        _, curr_sample, _ = get_file_name_pieces(curr_counts_fp)
        stripped_sample = strip_run_prefix(curr_sample, input_run_prefix)
        count_header, curr_counts_df = get_counts_df(curr_counts_fp, input_run_prefix)
        curr_counts_df.rename(columns={count_header:stripped_sample}, inplace=True)
        count_header = stripped_sample
        log2_series = make_log2_series(curr_counts_df[count_header], DEFAULT_PSEUDOCOUNT)
        
        title = " ".join([input_run_prefix, count_header, "with pseudocount", str(DEFAULT_PSEUDOCOUNT)])
        output_fp_prefix = build_multipart_fp(output_dir, [count_header, input_run_prefix])
        
        boxplot_fp = output_fp_prefix + "_" + boxplot_suffix
        show_and_save_boxplot(boxplot_fp, title, [count_header], log2_series)
        
        hist_fp = output_fp_prefix + "_" + "hist.png"
        show_and_save_histogram(hist_fp, title, log2_series)
        
        
def plot_combined_raw_counts(input_dir, input_run_prefix, combined_suffix, output_dir, output_run_prefix, boxplot_suffix):
    output_fp = build_multipart_fp(output_dir, [output_run_prefix, boxplot_suffix])
    combined_counts_fp = build_multipart_fp(input_dir, [input_run_prefix, combined_suffix])
    combined_counts_df = pandas.read_table(combined_counts_fp)
    samples_names = combined_counts_df.columns.values[1:]  # TODO: remove hardcode
    samples_data = []
    for curr_name in samples_names:
        log2_series = make_log2_series(combined_counts_df[curr_name], DEFAULT_PSEUDOCOUNT)
        samples_data.append(log2_series.tolist())
        
    title = " ".join([input_run_prefix, "all samples", "with pseudocount", str(DEFAULT_PSEUDOCOUNT)])
    show_and_save_boxplot(output_fp, title, samples_names, samples_data, 90)
