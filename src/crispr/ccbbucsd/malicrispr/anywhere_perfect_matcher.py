# standard libraries
import logging

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import rev_comp_canonical_dna_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def find_grna_names(grna_name_by_seq, fw_whole_seq, rv_whole_seq):
    grna_name_B = None
    grna_name_A = _find_grna_match_in_seq(grna_name_by_seq, fw_whole_seq)
    if grna_name_A is not None:
        rc_rv_whole_seq = rev_comp_canonical_dna_seq(rv_whole_seq)
        grna_name_B = _find_grna_match_in_seq(grna_name_by_seq, rc_rv_whole_seq)

    return grna_name_A, grna_name_B


def _find_grna_match_in_seq(names_and_seqs, larger_seq):
    found_names = []
    for potential_ref_name, potential_ref_seq in names_and_seqs:
        if larger_seq.find(potential_ref_seq) > -1:
            found_names.append(potential_ref_name)
            break

    if len(found_names) == 0:
        found_name = None
    else:
        found_name = found_names[-1]
        if len(found_names) > 1:
            logging.info("Matches for multiple names {0} found in sequence {1}".format(", ".join(found_names), larger_seq))
    # end if

    return found_name
