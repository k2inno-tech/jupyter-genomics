# standard libraries
import os

# third-party libraries
import nbformat
import nbparameterise
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor

# ccbb libraries
from ccbbucsd.utilities.files_and_paths import get_file_name_pieces, make_file_path

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def read_in_notebook(notebook_fp):
    with open(notebook_fp) as f:
        nb = nbformat.read(f, as_version=4)
    return nb


def set_parameters(nb, params_dict):
    orig_parameters = nbparameterise.extract_parameters(nb)
    params = nbparameterise.parameter_values(orig_parameters, **params_dict)
    new_nb = nbparameterise.replace_definitions(nb, params, execute=False)
    return new_nb


# modified from https://nbconvert.readthedocs.io/en/latest/execute_api.html
def execute_notebook(notebook_filename, notebook_filename_out, params_dict, run_path="", timeout=6000000):
    notebook_fp = os.path.join(run_path, notebook_filename)
    nb = read_in_notebook(notebook_fp)
    new_nb = set_parameters(nb, params_dict)
    ep = ExecutePreprocessor(timeout=timeout, kernel_name='python3')

    try:
        ep.preprocess(new_nb, {'metadata': {'path': run_path}})
    except:
        msg = 'Error executing the notebook "{0}".\n\n'.format(notebook_filename)
        msg = '{0}See notebook "{1}" for the traceback.'.format(msg, notebook_filename_out)
        print(msg)
        raise
    finally:
        with open(notebook_filename_out, mode='wt') as f:
            nbformat.write(new_nb, f)


def export_notebook_to_html(notebook_fp, output_dir):
    nb = read_in_notebook(notebook_fp)
    html_exporter = HTMLExporter()
    body, resources = html_exporter.from_notebook_node(nb)
    _, notebook_name, _ = get_file_name_pieces(notebook_fp)
    out_fp = make_file_path(output_dir, notebook_name, ".html")
    with open(out_fp, "w", encoding="utf8") as f:
        f.write(body)
