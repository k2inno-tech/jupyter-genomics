# standard libraries
import enum

# third-party libraries
import cutadapt.scripts.cutadapt

# ccbb libraries
from ccbbucsd.utilities.files_and_paths import get_file_name_pieces, make_file_path

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class TrimType(enum.Enum):
    FIVE = "5"
    THREE = "3"
    FIVE_THREE = "53"


def get_trimmed_suffix(trimtype):
    return "_trimmed{0}.fastq".format(trimtype.value)


def trim_linked_scaffold(output_dir, fastq_fp, scaffold_seq_5p, scaffold_seq_3p, quiet=True):
    args = ["-a", "{0}...{1}".format(scaffold_seq_5p,scaffold_seq_3p)]
    return _run_cutadapt(output_dir, fastq_fp, TrimType.FIVE_THREE, args, quiet)


def trim_global_scaffold(output_dir, fastq_fp, scaffold_seq_5p=None, scaffold_seq_3p=None, quiet=True):
    curr_fastq_fp = fastq_fp

    if scaffold_seq_5p is not None:
        curr_fastq_fp = _run_cutadapt_global(output_dir, curr_fastq_fp, scaffold_seq_5p, True, quiet)

    if scaffold_seq_3p is not None:
        curr_fastq_fp = _run_cutadapt_global(output_dir, curr_fastq_fp, scaffold_seq_3p, False, quiet)

    return curr_fastq_fp


def _run_cutadapt_global(output_dir, input_fastq_fp, seq_to_trim, is_5p, quiet):
    end_switch = "-g"
    end_name = TrimType.FIVE
    if not is_5p:
        end_switch = "-a"
        end_name = TrimType.THREE

    args = [end_switch, seq_to_trim]
    return _run_cutadapt(output_dir, input_fastq_fp, end_name, args, quiet)


def _run_cutadapt(output_dir, input_fastq_fp, trim_name, partial_args, quiet):
    _, input_base, _ = get_file_name_pieces(input_fastq_fp)
    output_fastq_fp = make_file_path(output_dir, input_base, get_trimmed_suffix(trim_name))
    args = [x for x in partial_args]
    if quiet:
        args.append("--quiet")
    args.extend(["-o", output_fastq_fp, input_fastq_fp])
    cutadapt.scripts.cutadapt.main(args)
    return output_fastq_fp
