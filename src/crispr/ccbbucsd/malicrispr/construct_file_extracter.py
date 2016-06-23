# standard libraries
import enum

# third-party libraries
import pandas

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import trim_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class ConstructHeaderKey(enum.Enum):
    PAST_INDEX = "past_index"
    CONSTRUCT_ID = "SequenceID"
    FINAL_SEQ = "FinalSequence"
    GENE_A_NAME = "Gene_A"
    GENE_A_CHR = "Gene_A_chr"
    GENE_A_POS = "Gene_A_pos"
    GENE_A_SEQ = "Gene_A_seq"
    GENE_B_NAME = "Gene_B"
    GENE_B_CHR = "Gene_B_chr"
    GENE_B_POS = "Gene_B_pos"
    GENE_B_SEQ = "Gene_B_seq"
    CONSTRUCT_A_NAME = "Construct_A_name"
    CONSTRUCT_B_NAME = "Construct_B_name"


def get_construct_separator():
    return "__"


def extract_construct_and_grna_info(constructs_fp):
    construct_table = _read_in_construct_table(constructs_fp, rows_to_skip=1)
    seq_name_sets = _extract_grnas_from_construct_table(construct_table)
    grna_name_seq_pairs = _format_and_check_grnas_input(seq_name_sets)
    construct_names = construct_table[ConstructHeaderKey.CONSTRUCT_ID.value].unique().tolist()
    return construct_names, grna_name_seq_pairs


def _read_in_construct_table(constructs_fp, rows_to_skip=1):
    col_names = [x.value for x in ConstructHeaderKey]
    return pandas.read_table(constructs_fp, skiprows=rows_to_skip, names=col_names)


def _extract_grnas_from_construct_table(construct_table):
    grna_name_key = "grna_name"
    grna_seq_key = "grna_seq"

    # split the construct id in each row into two pieces--the two construct names--and put them into new columns
    construct_table[ConstructHeaderKey.CONSTRUCT_A_NAME.value], \
    construct_table[ConstructHeaderKey.CONSTRUCT_B_NAME.value] = \
        zip(*construct_table[ConstructHeaderKey.CONSTRUCT_ID.value].str.split(get_construct_separator()).tolist())

    # get the gene/sequence pairs for each of the two genes and concatenate them
    gene_a_pairs = _extract_grna_name_and_seq_df(construct_table, ConstructHeaderKey.CONSTRUCT_A_NAME.value,
                                                 ConstructHeaderKey.GENE_A_SEQ.value, grna_name_key, grna_seq_key)
    gene_b_pairs = _extract_grna_name_and_seq_df(construct_table, ConstructHeaderKey.CONSTRUCT_B_NAME.value,
                                                 ConstructHeaderKey.GENE_B_SEQ.value, grna_name_key, grna_seq_key)
    gene_pairs = pandas.concat([gene_a_pairs, gene_b_pairs])

    # extract only the unique pairs
    grna_seq_name_groups = gene_pairs.groupby([grna_seq_key, grna_name_key]).groups
    return [x for x in grna_seq_name_groups]


def _extract_grna_name_and_seq_df(construct_table, name_key, seq_key, grna_name_key, grna_seq_key):
    name_and_seq_df = construct_table[[name_key, seq_key]]
    name_and_seq_df.rename(columns={name_key: grna_name_key, seq_key: grna_seq_key}, inplace=True)
    return name_and_seq_df


def trim_grnas(grnas_name_and_seq_list, retain_len):
    result = []
    for name_seq_tuple in grnas_name_and_seq_list:
        grna_name = name_seq_tuple[0]
        full_seq = name_seq_tuple[1]
        trimmed_seq = trim_seq(full_seq, retain_len, False)  # False = do not retain from 5p end but from 3p end
        result.append((grna_name, trimmed_seq))
    return result


def _read_in_grnas(grnas_fp):
    list_from_file = _read_grnas_input(grnas_fp)
    return _format_and_check_grnas_input(list_from_file)


def _read_grnas_input(grnas_fp, comment_prefix="#", delimiter="\t"):
    result = []
    with open(grnas_fp, 'r') as file_handle:
        for line in file_handle:
            if not line.startswith(comment_prefix):
                pieces = line.strip().split(delimiter)
                result.append(pieces)
    return result


def _format_and_check_grnas_input(grnas_seq_and_name_list):
    expected_num_pieces = 2
    seqs_by_names = {}
    names_by_seqs = {}
    result = []

    for curr_set in grnas_seq_and_name_list:
        if len(curr_set) != expected_num_pieces:
            raise ValueError(
                "input '{0}' has {1} pieces instead of the expected {2}".format(
                    curr_set, len(curr_set), expected_num_pieces
                ))
        curr_seq = curr_set[0]
        curr_name = curr_set[1]

        if curr_seq in names_by_seqs:
            raise ValueError(
                "sequence '{0}' associated with name '{1}' but was already associated with name '{2}'".format(
                    curr_seq, curr_name, names_by_seqs[curr_seq]
                ))

        if curr_name in seqs_by_names:
            raise ValueError(
                "name '{0}' associated with sequence '{1}' but was already associated with sequence '{2}'".format(
                    curr_name, curr_seq, seqs_by_names[curr_name]
                ))

        names_by_seqs[curr_seq] = curr_name
        seqs_by_names[curr_name] = curr_seq
        result.append((curr_name, curr_seq))
    # next pair in

    return result
