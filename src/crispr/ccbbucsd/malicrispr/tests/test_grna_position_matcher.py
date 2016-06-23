# standard libraries
import unittest

# project-specific libraries
from ccbbucsd.malicrispr.grna_position_matcher import GrnaPositionMatcher

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestGrnaPositionMatcher(unittest.TestCase):
    # region _id_sequence_match tests
    def test__id_sequence_match_perfect(self):
        expected_name = "test_grna_2"
        test_seq = "AAAT"
        names_and_seqs = [("test_grna_1", "ACCG"),
                          (expected_name, test_seq)]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)

        output = matcher._id_sequence_match(1, test_seq)
        self.assertEqual(expected_name, output)

    def test__id_sequence_match_acceptable_mismatch(self):
        expected_name = "test_grna_2"
        names_and_seqs = [("test_grna_1", "ACCG"),
                          (expected_name, "AAAT")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)

        output = matcher._id_sequence_match(1, "AACT")
        self.assertEqual(expected_name, output)

    def test__id_sequence_match_unacceptable_mismatch(self):
        expected_name = "test_grna_2"
        names_and_seqs = [("test_grna_1", "ACCG"),
                          (expected_name, "AAAT")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)

        output = matcher._id_sequence_match(1, "TTCA")
        self.assertIsNone(output)

    # endregion

    # region _id_pair_matches tests
    def test__id_pair_matches_both(self):
        name1 = "test_grna_1"
        name2 = "test_grna_2"
        names_and_seqs = [(name1, "ACCG"),
                          (name2, "AAAT")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)
        grnaA, grnaB = matcher._id_pair_matches("AAAA", "ACCC")
        self.assertEqual(name2, grnaA)
        self.assertEqual(name1, grnaB)

    def test__id_pair_matches_a_but_not_b(self):
        name1 = "test_grna_1"
        name2 = "test_grna_2"
        names_and_seqs = [(name1, "ACCG"),
                          (name2, "AAAT")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)
        grnaA, grnaB = matcher._id_pair_matches("AAAA", "AGGG")
        self.assertEqual(name2, grnaA)
        self.assertIsNone(grnaB)

    def test__id_pair_matches_neither(self):
        name1 = "test_grna_1"
        name2 = "test_grna_2"
        names_and_seqs = [(name1, "ACCG"),
                          (name2, "AAAT")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)
        grnaA, grnaB = matcher._id_pair_matches("TTAA", "AGGG")
        self.assertIsNone(grnaA)
        self.assertIsNone(grnaB)

    # endregion

    # region _generate_seqs_to_check tests
    def test__generate_seqs_to_check(self):
        fw_output, rv_output = GrnaPositionMatcher._generate_seqs_to_check("ACCG", "AAAT")
        self.assertEqual("ACCG", fw_output)
        self.assertEqual("ATTT", rv_output)

    # endregion

    # region find_fw_and_rv_read_matches tests
    def test_find_fw_and_rv_read_matches(self):
        name2 = "test_grna_2"
        name3 = "something_else"
        names_and_seqs = [("test_grna_1", "ACCG"),
                          (name2, "AAAT"),
                          (name3, "GGAG")]
        matcher = GrnaPositionMatcher(names_and_seqs, 4)
        grnaA, grnaB = matcher.find_fw_and_rv_read_matches("AACT", "CTCA")
        self.assertEqual(grnaA, name2)
        self.assertEqual(grnaB, name3)

    # end region