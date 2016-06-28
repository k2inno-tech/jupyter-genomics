# standard libraries
import enum
import os

# ccbb libraries
from ccbbucsd.utilities.notebook_pipeliner import DATASET_NAME_KEY, ALG_NAME_KEY, execute_run

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class PipelineSteps(enum.Enum):
    SCAFFOLD_TRIMMING = 0
    TRIMMED_READ_FILTERING = 1
    CONSTRUCT_COUNTING = 2
    COUNT_COMBINATION = 3
    COUNT_PLOTTING = 4


def main():
    # *********************************************************
    # you'll probably change these on every run, unless you're reprocessing an existing dataset for some reason
    fastq_set_name = "20160627_HeLa_A549_CV4"
    human_readable_name = "20160627_HeLa_A549_CV4"

    # *********************************************************
    # you'll probably change these lines only when changing the construct library from which you're analyzing data.
    # Note that the file named here is expected to be in tab-delimited text format with no quotes around text values,
    # and that the first line will be skipped (as a header row).

    # option 1: info for CV4
    spacers_file_name = "CV4_2spacers.txt"  # for CV4
    full_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
    full_5p_r2 = "CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"
    full_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGC"
    full_3p_r2 = "CAAACAAGGCTTTTCTCCAAGG"

    # option 2: info for MV4
    # spacers_file_name = "Metabolism_dual_spacers.txt"
    # full_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
    # full_5p_r2 = "CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"
    # full_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG"
    # full_3p_r2 = "CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA"

    # *********************************************************
    # you'll probably decide on these values before beginning screen analysis, set them once, and then leave them alone
    # unless you change your analysis approach and/or hardware
    len_of_seq_to_match = 19
    num_allowed_mismatches = 1
    num_processors = 3
    main_dir = "/home/ec2-user/jupyter-genomics"

    # *********************************************************
    # you'll probably never change these, assuming you use the suggested file structure
    notebooks_dir = os.path.join(main_dir, "notebooks/crispr")
    data_dir = "/data"
    raw_dir = os.path.join(data_dir, "raw")
    library_def_dir = os.path.join(main_dir, "library_definitions")
    interim_fastq_set_dir = os.path.join(data_dir, "interim", fastq_set_name)
    processed_dir = os.path.join(data_dir, "processed")

    # *********************************************************
    # DON'T change these unless you *really* know what you're doing AND either you're restarting a run from part-way
    # through or you're refactoring the notebook parameters
    steps_to_run = [PipelineSteps.SCAFFOLD_TRIMMING,
                    PipelineSteps.TRIMMED_READ_FILTERING,
                    PipelineSteps.CONSTRUCT_COUNTING,
                    PipelineSteps.COUNT_COMBINATION,
                    PipelineSteps.COUNT_PLOTTING]

    shared_params = {'g_code_location': os.path.join(main_dir, "src/crispr"),
                     'g_timestamp': "{timestamp}",
                     DATASET_NAME_KEY: human_readable_name,
                     'g_num_processors': num_processors,
                     'g_fastqs_dir': os.path.join(raw_dir, fastq_set_name),
                     'g_trimmed_fastqs_dir': interim_fastq_set_dir,
                     'g_full_5p_r1': full_5p_r1,
                     'g_full_5p_r2': full_5p_r2,
                     'g_full_3p_r1': full_3p_r1,
                     'g_full_3p_r2': full_3p_r2,
                     'g_filtered_fastqs_dir': interim_fastq_set_dir,
                     ALG_NAME_KEY: "{0}mer_{1}mm_py".format(len_of_seq_to_match, num_allowed_mismatches),
                     'g_len_of_seq_to_match': len_of_seq_to_match,
                     'g_num_allowed_mismatches': num_allowed_mismatches,
                     'g_constructs_fp': os.path.join(library_def_dir, spacers_file_name),
                     'g_fastq_counts_dir': interim_fastq_set_dir,
                     'g_fastq_counts_run_prefix': human_readable_name + "_{timestamp}",
                     'g_collapsed_counts_dir': "{run_dir}",
                     'g_collapsed_counts_run_prefix': "{run_prefix}",
                     'g_combined_counts_dir': "{run_dir}",
                     'g_combined_counts_run_prefix': "{run_prefix}",
                     'g_plots_dir': "{run_dir}",
                     'g_plots_run_prefix': "{run_prefix}"
                     }

    # *********************************************************
    # DON'T CHANGE ANYTHING BELOW THIS LINE UNLESS YOU ARE AMANDA BIRMINGHAM.  And even then think twice :)
    pipeline_steps = {PipelineSteps.SCAFFOLD_TRIMMING: [notebooks_dir, "Dual CRISPR 1-Construct Scaffold Trimming.ipynb"],
                      PipelineSteps.TRIMMED_READ_FILTERING: [notebooks_dir, "Dual CRISPR 2-Constuct Filter.ipynb"],
                      PipelineSteps.CONSTRUCT_COUNTING: [notebooks_dir, "Dual CRISPR 3-Construct Counting.ipynb"],
                      PipelineSteps.COUNT_COMBINATION: [notebooks_dir, "Dual CRISPR 4-Count Combination.ipynb"],
                      PipelineSteps.COUNT_PLOTTING: [notebooks_dir, "Dual CRISPR 5-Count Plots.ipynb"]}

    execute_run(pipeline_steps, shared_params, steps_to_run, processed_dir)


if __name__ == '__main__':
    main()
