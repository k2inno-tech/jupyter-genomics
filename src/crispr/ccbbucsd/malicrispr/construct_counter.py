"""This module counts almost-perfect matches of small sequences within forward and reverse fastq sequence pairs."""

# standard libraries
import csv
import datetime
import logging

# ccbb libraries
from ccbbucsd.utilities.basic_fastq import FastqHandler, paired_fastq_generator

# project-specific libraries
from ccbbucsd.malicrispr.construct_file_extracter import get_construct_separator

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def get_counts_file_suffix():
    return "counts.txt"


def get_counter_from_names(names_to_count):
    return {x: 0 for x in names_to_count}


def generate_construct_counts(grna_matcher, construct_names, output_fp, fw_fastq_fp, rv_fastq_fp):
    counts_info_tuple = _match_and_count_constructs_from_files(grna_matcher, construct_names, fw_fastq_fp, rv_fastq_fp)
    counts_by_construct = counts_info_tuple[0]
    counts_by_type = counts_info_tuple[1]
    _write_counts(counts_by_construct, counts_by_type, output_fp)


def _match_and_count_constructs_from_files(grna_matcher, construct_names, fw_fastq_fp, rv_fastq_fp):
    construct_counts = get_counter_from_names(construct_names)
    fw_fastq_handler = FastqHandler(fw_fastq_fp)
    rv_fastq_handler = FastqHandler(rv_fastq_fp)
    return _match_and_count_constructs(grna_matcher, construct_counts, fw_fastq_handler, rv_fastq_handler)


def _match_and_count_constructs(grna_matcher, construct_counts, fw_fastq_handler, rv_fastq_handler):
    summary_counts = {"num_pairs": 0, "num_pairs_unrecognized": 0, "num_constructs_found": 0,
                      "num_constructs_unrecognized": 0, "num_constructs_recognized": 0}

    paired_fastq_seqs = paired_fastq_generator(fw_fastq_handler, rv_fastq_handler)
    for curr_pair_seqs in paired_fastq_seqs:
        summary_counts["num_pairs"] += 1
        _report_progress(summary_counts["num_pairs"])

        grna_name_A, grna_name_B = grna_matcher.find_fw_and_rv_read_matches(*curr_pair_seqs)
        if grna_name_A is not None and grna_name_B is not None:
            construct_name = _generate_construct_name(grna_name_A, grna_name_B)
            summary_counts["num_constructs_found"] += 1

            if construct_name in construct_counts:
                summary_counts["num_constructs_recognized"] += 1
                construct_counts[construct_name] += 1
            else:
                summary_counts["num_constructs_unrecognized"] += 1
                logging.debug("Unrecognized construct name: {0}".format(construct_name))
        else:
            summary_counts["num_pairs_unrecognized"] += 1
            logging.debug("Unrecognized sequence: {0},{1}".format(*curr_pair_seqs))

    return construct_counts, summary_counts


def _report_progress(num_fastq_pairs):
    if num_fastq_pairs % 100000 == 0:
        logging.info("On fastq pair number {0} at {1}".format(num_fastq_pairs, datetime.datetime.now()))


def _generate_construct_name(grna_name_A, grna_name_B):
    return "{0}{1}{2}".format(grna_name_A, get_construct_separator(), grna_name_B)


def _write_counts(counts_by_construct, counts_by_type, output_fp):
    construct_names = list(counts_by_construct.keys())
    construct_names.sort()

    with open(output_fp, 'w') as file_handle:
        summary_pieces = []
        for curr_key, curr_value in counts_by_type.items():
            summary_pieces.append("{0}:{1}".format(curr_key, curr_value))
        summary_comment = ",".join(summary_pieces)
        summary_comment = "# " + summary_comment
        header = ["construct_id", "counts"]
        writer = csv.writer(file_handle, delimiter="\t")

        writer.writerow([summary_comment])
        writer.writerow(header)
        for curr_name in construct_names:
            row = [curr_name, counts_by_construct[curr_name]]
            writer.writerow(row)