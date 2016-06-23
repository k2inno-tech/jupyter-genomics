# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import rev_comp_canonical_dna_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class GrnaPositionMatcher:
    @staticmethod
    def _generate_seqs_to_check(fw_whole_seq, rv_whole_seq):
        rc_whole_rv_seq = rev_comp_canonical_dna_seq(rv_whole_seq)
        return fw_whole_seq, rc_whole_rv_seq

    def __init__(self, grna_names_and_seqs, expected_len, num_allowed_fw_mismatches, num_allowed_rv_mismatches):
        self._grna_names_and_seqs = grna_names_and_seqs
        self._num_allowed_fw_mismatches = num_allowed_fw_mismatches
        self._num_allowed_rv_mismatches = num_allowed_rv_mismatches
        self._seq_len = expected_len

    @property
    def num_allowed_fw_mismatches(self):
        return self._num_allowed_fw_mismatches

    @property
    def num_allowed_rv_mismatches(self):
        return self._num_allowed_rv_mismatches

    def find_fw_and_rv_read_matches(self, fw_whole_seq, rv_whole_seq):
        fw_construct_window, rc_rv_construct_window = self._generate_seqs_to_check(fw_whole_seq, rv_whole_seq)
        return self._id_pair_matches(fw_construct_window, rc_rv_construct_window)

    def _id_pair_matches(self, input1_seq, input2_seq):
        input2_match_name = None
        input1_match_name = self._id_sequence_match(self.num_allowed_fw_mismatches, input1_seq)
        if input1_match_name is not None:
            input2_match_name = self._id_sequence_match(self.num_allowed_rv_mismatches, input2_seq)

        return input1_match_name, input2_match_name

    def _id_sequence_match(self, num_allowed_mismatches, input_seq):
        found_name = None

        # check for perfect matches
        for potential_name, potential_reference in self._grna_names_and_seqs:
            if input_seq == potential_reference:
                found_name = potential_name
                break

        # if no perfect matches found, check for mismatches
        if found_name is None:
            min_found_mismatches = num_allowed_mismatches + 1  # nothing checked yet so num mismatches is maximum
            found_name = None
            for potential_name, potential_reference in self._grna_names_and_seqs:

                # this way is not elegant, but it IS fast :) .... >3x faster than XOR approach
                num_mismatches = 0
                for x in range(0, self._seq_len):
                    if input_seq[x] != potential_reference[x]:  # NB assumption that both are the same fixed length!
                        num_mismatches += 1
                        if num_mismatches > num_allowed_mismatches:
                            break

                if num_mismatches < min_found_mismatches:
                    min_found_mismatches = num_mismatches
                    if num_mismatches <= num_allowed_mismatches:
                        found_name = potential_name

        return found_name
