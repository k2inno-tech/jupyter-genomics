# standard libraries
import logging

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import trim_seq
from ccbbucsd.utilities.basic_fastq import FastqHandler, paired_fastq_generator
from ccbbucsd.utilities.files_and_paths import transform_path

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def get_filtered_file_suffix():
    return "_len_filtered.fastq"


def filter_pair_by_len(min_len, max_len, retain_len, output_dir, fw_fastq_fp, rv_fastq_fp):
    fw_fastq_handler = FastqHandler(fw_fastq_fp)
    rv_fastq_handler = FastqHandler(rv_fastq_fp)
    fw_out_handle, rv_out_handle = _open_output_file_pair(fw_fastq_fp, rv_fastq_fp, output_dir)
    counters = {"num_pairs": 0, "num_pairs_passing": 0}

    filtered_fastq_records = _filtered_fastq_generator(fw_fastq_handler, rv_fastq_handler, min_len, max_len, retain_len,
                                                       counters)
    for fw_record, rv_record in filtered_fastq_records:
        fw_out_handle.writelines(fw_record.lines)
        rv_out_handle.writelines(rv_record.lines)

    fw_out_handle.close()
    rv_out_handle.close()
    return _summarize_counts(counters)


def _filtered_fastq_generator(fw_fastq_handler, rv_fastq_handler, min_len, max_len, retain_len, counters):
    paired_fastq_records = paired_fastq_generator(fw_fastq_handler, rv_fastq_handler, True)
    for curr_pair_fastq_records in paired_fastq_records:
        counters["num_pairs"] += 1
        _report_progress(counters["num_pairs"])

        fw_record = curr_pair_fastq_records[0]
        fw_passing_seq = _check_and_trim_seq(_get_upper_seq(fw_record), min_len, max_len, retain_len, False)
        if fw_passing_seq is not None:
            rv_record = curr_pair_fastq_records[1]
            rv_passing_seq = _check_and_trim_seq(_get_upper_seq(rv_record), min_len, max_len, retain_len, True)
            if rv_passing_seq is not None:
                counters["num_pairs_passing"] += 1
                fw_record.sequence = fw_passing_seq
                fw_record.quality = trim_seq(fw_record.quality, retain_len, False)
                rv_record.sequence = rv_passing_seq
                rv_record.quality = trim_seq(rv_record.quality, retain_len, True)
                yield fw_record, rv_record


def _open_output_file_pair(fw_fastq_fp, rv_fastq_fp, output_dir):
    fw_fp = transform_path(fw_fastq_fp, output_dir, get_filtered_file_suffix())
    rv_fp = transform_path(rv_fastq_fp, output_dir, get_filtered_file_suffix())
    fw_handle = open(fw_fp, 'w')
    rv_handle = open(rv_fp, 'w')
    return fw_handle, rv_handle


def _report_progress(num_fastq_pairs):
    if num_fastq_pairs % 100000 == 0:
        logging.debug("On fastq pair number {0}".format(num_fastq_pairs))


def _get_upper_seq(fastq_record):
    return fastq_record.sequence.upper()


def _check_and_trim_seq(input_seq, min_len, max_len, retain_len, retain_5p_end):
    result = None
    seq_len = len(input_seq)
    if seq_len >= min_len and seq_len <= max_len:
        result = trim_seq(input_seq, retain_len, retain_5p_end)
    return result


# def _write_fasta_record(file_handle, label, sequence):
#     lines = [">" + label, sequence]
#     file_handle.writelines(lines)


def _summarize_counts(counts_by_type):
    summary_pieces = []
    sorted_keys = sorted(counts_by_type.keys())  # sort to ensure deterministic output ordering
    for curr_key in sorted_keys:
        curr_value = counts_by_type[curr_key]
        summary_pieces.append("{0}:{1}".format(curr_key, curr_value))
    result = ",".join(summary_pieces)
    return result
