# standard libraries
import os

# ccbb libraries
from ccbbucsd.utilities.files_and_paths import get_wild_path, group_files

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def rev_comp_canonical_dna_seq(dna_seq):
    complement_dict = {"A": "T",
                       "C": "G",
                       "G": "C",
                       "T": "A",
                       "a": "t",
                       "c": "g",
                       "g": "c",
                       "t": "a"}
    reversed_seq = dna_seq[::-1]
    rev_complement_chars = [complement_dict[x] if x in complement_dict else x for x in reversed_seq]
    return "".join(rev_complement_chars)


def expand_possible_mismatches(perfect_seq, position_alphabet, include_perfect=False):
    result = []
    if include_perfect:
        result.append(perfect_seq)

    perfect_seq_chars = perfect_seq.split()
    for curr_position in range(0, len(perfect_seq_chars)):
        perfect_seq_char = perfect_seq_chars[curr_position]
        for possible_char in position_alphabet:
            if possible_char != perfect_seq_char:
                mismatch_seq_chars = list(perfect_seq_chars)
                mismatch_seq_chars[curr_position] = possible_char
                result.append("".join(mismatch_seq_chars))

    return result


def pair_hiseq_read_files(fastq_filepaths):
    paired_fastqs_by_base = group_files(fastq_filepaths, "_R\d_", "_")

    num_expected = 2
    failure_msgs = []
    for curr_key in paired_fastqs_by_base:
        num_paths = len(paired_fastqs_by_base[curr_key])

        if num_paths != num_expected:
            failure_msgs.append(
                "{0} has {1} read files instead of {2}".format(curr_key, num_paths, num_expected))

    if len(failure_msgs) == 0:
        failure_msgs = None
    else:
        failure_msgs = "\n".join(failure_msgs)

    return paired_fastqs_by_base, failure_msgs


def gunzip_fastqs(directory):
    # gunzip the gzipped files; do this from shell because doing through
    # python gzip module is slow
    gz_wildpath = get_wild_path(directory, '*.fastq.gz')
    gunzip_cmd = "gunzip -k " + gz_wildpath
    os.system(gunzip_cmd)


def trim_seq(input_seq, retain_len, retain_5p_end):
    if len(input_seq) < retain_len:
        raise ValueError(
            "input sequence {0} has length {1}, shorter than retain length {2}".format(input_seq, len(input_seq),
                                                                                       retain_len))

    if retain_5p_end:
        return input_seq[:retain_len]
    else:
        return input_seq[-retain_len:]


