# standard libraries
import os
import re

# ccbb libraries
from ccbbucsd.utilities.analysis_run_prefixes import get_timestamp, get_run_prefix
from ccbbucsd.utilities.files_and_paths import verify_or_make_dir
from ccbbucsd.utilities.notebook_runner import execute_notebook, export_notebook_to_html

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


DATASET_NAME_KEY = "g_dataset_name"
ALG_NAME_KEY = "g_count_alg_name"


def execute_run(possible_actions_dict, run_params, ordered_run_steps, parent_dir, run_folder=None):
    timestamp = get_timestamp()
    run_prefix = _generate_run_prefix(run_params, timestamp)
    if run_folder is None:
        run_dir = _make_run_dir(parent_dir, timestamp)
    else:
        run_dir = os.path.join(parent_dir, run_folder)
    methods_dir = _create_run_and_methods_dirs(run_dir)

    for run_action in ordered_run_steps:
        step_details = possible_actions_dict[run_action]
        run_and_output_notebook(step_details, run_params, timestamp, run_prefix, run_dir, methods_dir)


def run_and_output_notebook(step_settings, params_dict, timestamp, run_prefix, run_dir, methods_dir):
    run_path = step_settings[0]
    base_notebook_filename = step_settings[1]

    formatted_params_dict = _format_parameters(run_dir, timestamp, run_prefix, params_dict)
    notebook_out_fp = _get_output_fp(base_notebook_filename, timestamp, methods_dir, ".ipynb")
    execute_notebook(base_notebook_filename, notebook_out_fp, formatted_params_dict, run_path)
    export_notebook_to_html(notebook_out_fp, methods_dir)


def _generate_run_prefix(run_params, timestamp):
    dataset_name = run_params[DATASET_NAME_KEY]
    alg_name = run_params[ALG_NAME_KEY]
    return get_run_prefix(dataset_name, alg_name, timestamp)


def _make_run_dir(parent_path, timestamp):
    return os.path.join(parent_path, "run{0}".format(timestamp))


def _create_run_and_methods_dirs(run_dir_name):
    methods_dir_name = os.path.join(run_dir_name, _get_methods_folder_name())
    verify_or_make_dir(run_dir_name)
    verify_or_make_dir(methods_dir_name)
    return methods_dir_name


def _get_methods_folder_name():
    return "methods"


def _mangle_notebook_name(timestamp, notebook_filename):
    delimiter = "_"
    name_base, _ = os.path.splitext(notebook_filename)
    lower_base = name_base.lower()
    delimited_notebook_base = re.sub("\s+", delimiter, lower_base)
    new_base = "{0}{1}{2}".format(timestamp, delimiter, delimited_notebook_base)
    return new_base


def _get_output_fp(notebook_fp, timestamp, methods_dir, output_ext):
    _, notebook_name = os.path.split(notebook_fp)
    new_base = _mangle_notebook_name(timestamp, notebook_name)
    new_fp = os.path.join(methods_dir, new_base + output_ext)
    return new_fp


def _format_parameters(output_dir, timestamp, run_prefix, params_dict):
    result = {}
    for curr_param_key, curr_param_item in params_dict.items():
        if hasattr(curr_param_item, 'format'):  # if this item has a format method
            final_item = curr_param_item.format(run_dir=output_dir, timestamp=timestamp, run_prefix=run_prefix)
        else:
            final_item = curr_param_item
        result[curr_param_key] = final_item
    return result


