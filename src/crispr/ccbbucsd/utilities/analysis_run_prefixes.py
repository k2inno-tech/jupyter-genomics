# standard libraries
import datetime

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_timestamp():
    return datetime.datetime.now().strftime('%Y%m%d%H%M%S')


def join_name_pieces(name_pieces):
    return "_".join(name_pieces)    


def get_run_prefix(dataset_name, alg_name, timestamp):
    pieces = [dataset_name, alg_name, timestamp]
    non_empty_pieces = [x for x in pieces if (x is not None) and (x.strip())]
    return join_name_pieces(non_empty_pieces)


def strip_run_prefix(string_to_strip, run_prefix):
    run_prefix_with_separator = join_name_pieces([run_prefix, ""])  # empty string to force ending separator
    result = string_to_strip.replace(run_prefix_with_separator, "")  # remove w/ separator
    result = result.replace(run_prefix, "")  # if exists w/o separator, remove that too
    return result


def check_or_set(input_val, output_val):
    return input_val if input_val else output_val


def describe_var_list(input_var_name_list):
    description_list =  ["{0}: {1}\n".format(name, eval(name)) for name in input_var_name_list]
    return "".join(description_list)