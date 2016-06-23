# ccbb libraries
from ccbbucsd.utilities.analysis_run_prefixes import strip_run_prefix
from ccbbucsd.utilities.files_and_paths import build_multipart_fp, group_files, get_filepaths_by_prefix_and_suffix

# project-specific libraries
from ccbbucsd.malicrispr.count_files_and_dataframes import get_counts_df

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_collapsed_counts_file_suffix():
    return "collapsed.txt"


def get_combined_counts_file_suffix():
    return "counts_combined.txt"


def group_lane_and_set_files(filepaths):
    # NB: this regex assumes read designator has *already* been removed
    # and replaced with _ as done by group_read_pairs
    return group_files(filepaths, "_L\d\d\d_\d\d\d", "")


def combine_count_files(counts_fp_for_dataset, run_prefix):
    combined_df = None
    
    for curr_counts_fp in counts_fp_for_dataset:
        count_header, curr_counts_df = get_counts_df(curr_counts_fp, run_prefix)
        
        if combined_df is None:
            combined_df = curr_counts_df
        else:
            combined_df[count_header] = curr_counts_df[count_header]
    
    return combined_df


def write_collapsed_count_files(input_dir, output_dir, curr_run_prefix, counts_run_prefix, counts_suffix, counts_collapsed_file_suffix):
    counts_fps_for_dataset = get_filepaths_by_prefix_and_suffix(input_dir, counts_run_prefix, counts_suffix)
    fps_by_sample = group_lane_and_set_files(counts_fps_for_dataset)
    
    for curr_sample, curr_fps in fps_by_sample.items():
        stripped_sample = strip_run_prefix(curr_sample, counts_run_prefix)
        output_fp = build_multipart_fp(output_dir, [curr_run_prefix, stripped_sample, counts_collapsed_file_suffix]) 
        combined_df = None        
        
        for curr_fp in curr_fps:
            count_header, curr_counts_df = get_counts_df(curr_fp, counts_run_prefix)
        
            if combined_df is None:
                combined_df = curr_counts_df
                combined_df.rename(columns = {count_header:stripped_sample}, inplace = True) 
            else:
                combined_df[stripped_sample] = combined_df[stripped_sample] + curr_counts_df[count_header]
    
        combined_df.to_csv(output_fp, sep="\t", index=False)    


def write_combined_count_file(input_dir, output_dir, curr_run_prefix, counts_run_prefix, counts_suffix, combined_suffix):
    output_fp = build_multipart_fp(output_dir, [curr_run_prefix, combined_suffix])
    counts_fps_for_run = get_filepaths_by_prefix_and_suffix(input_dir, counts_run_prefix, counts_suffix)
    combined_df = combine_count_files(counts_fps_for_run, curr_run_prefix)
    combined_df.to_csv(output_fp, sep="\t", index=False)