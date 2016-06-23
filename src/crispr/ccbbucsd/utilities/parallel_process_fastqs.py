# standard libraries
import datetime
import logging
import multiprocessing
import timeit
import traceback

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import pair_hiseq_read_files
from ccbbucsd.utilities.files_and_paths import get_filepaths_from_wildcard

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_elapsed_time_to_now(start_time, process_name=None):
    end_time = timeit.default_timer()
    elapsed_seconds = end_time - start_time
    m, s = divmod(elapsed_seconds, 60)
    h, m = divmod(m, 60)
    result = "elapsed time: %d:%02d:%02d" % (h, m, s)
    if process_name is not None:
        result = "{0} ".format(process_name) + result
    return result


def time_function(process_name, func_name, pass_process_name_to_func, *func_args):
    logging.info("Starting {0} at {1}".format(process_name, datetime.datetime.now()))
    start_time = timeit.default_timer()

    if pass_process_name_to_func:
        try:
            func_result = func_name(process_name, *func_args)
        except:
            func_result = traceback.format_exc()
    else:
        try:
            func_result = func_name(*func_args)
        except:
            func_result = traceback.format_exc()

    logging.info(get_elapsed_time_to_now(start_time, process_name))
    return process_name, func_result


def parallel_process_paired_reads(fastq_dir, file_suffix, num_processes, func_for_one_pair, func_fixed_inputs_list,
                                  pass_process_name_to_func=False):
    logging.info("Starting parallel processing at {0}".format(datetime.datetime.now()))
    start_time = timeit.default_timer()

    results = []
    fastq_filepaths = get_filepaths_from_wildcard(fastq_dir, file_suffix)
    paired_fastqs_by_base, failure_msgs = pair_hiseq_read_files(fastq_filepaths)

    if failure_msgs is not None:
        logging.info(failure_msgs)
    else:
        process_arguments = []
        sorted_bases = sorted(paired_fastqs_by_base.keys())
        for curr_base in sorted_bases:
            fp_list = paired_fastqs_by_base[curr_base]
            curr_args_list = [curr_base, func_for_one_pair, pass_process_name_to_func]
            curr_args_list.extend(func_fixed_inputs_list)
            curr_args_list.extend(fp_list)
            process_arguments.append(tuple(curr_args_list))

        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.starmap(time_function, process_arguments)

    logging.info(get_elapsed_time_to_now(start_time, "parallel processing"))
    return results


def concatenate_parallel_results(results_tuples):
    results_lines = ["{0}: {1}\n".format(x[0], x[1] if x[1] is not None else "finished")
                     for x in results_tuples]
    return "".join(results_lines)
