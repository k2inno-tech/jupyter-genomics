# standard libraries
import unittest

# project-specific libraries
from ccbbucsd.malicrispr.count_filterer import _check_and_trim_seq, _filtered_fastq_generator, _summarize_counts, \
    _trim_seq
from ccbbucsd.utilities.basic_fastq import FastqHandler

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region _trim_seq tests
    def test__trim_seq_long(self):
        input_seq = "ACGT"
        retain_len = 3

        # trim from 5p end
        output_5p = _trim_seq(input_seq, retain_len, False)
        self.assertEqual("CGT", output_5p)

        # trim from 3p end
        output_3p = _trim_seq(input_seq, retain_len, True)
        self.assertEqual("ACG", output_3p)

    def test__trim_seq_exact(self):
        input_seq = "ACGT"
        retain_len = 4

        # trim from 5p end
        output_5p = _trim_seq(input_seq, retain_len, False)
        self.assertEqual(input_seq, output_5p)

        # trim from 3p end
        output_3p = _trim_seq(input_seq, retain_len, True)
        self.assertEqual(input_seq, output_3p)

    def test__trim_seq_short(self):
        input_seq = "ACGT"
        retain_len = 5

        # trim from 5p end
        with self.assertRaises(ValueError):
            _trim_seq(input_seq, retain_len, False)

        # trim from 3p end
        with self.assertRaises(ValueError):
            _trim_seq(input_seq, retain_len, True)

    # endregion

    # region _check_and_trim_seq tests
    def test__check_and_trim_seq_short(self):
        input_seq = "ACGT"
        min_len = 5
        max_len = 8
        retain_len = 7

        # trim from 5p end
        output_5p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, False)
        self.assertIsNone(output_5p)

        # trim from 3p end
        output_3p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, True)
        self.assertIsNone(output_3p)

    def test__check_and_trim_seq_acceptable(self):
        input_seq = "ACGT"
        min_len = 3
        max_len = 4
        retain_len = 3

        # trim from 5p end
        output_5p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, False)
        self.assertEqual("CGT", output_5p)

        # trim from 3p end
        output_3p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, True)
        self.assertEqual("ACG", output_3p)

    def test__check_and_trim_seq_long(self):
        input_seq = "ACGT"
        min_len = 1
        max_len = 3
        retain_len = 2

        # trim from 5p end
        output_5p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, False)
        self.assertIsNone(output_5p)

        # trim from 3p end
        output_3p = _check_and_trim_seq(input_seq, min_len, max_len, retain_len, True)
        self.assertIsNone(output_3p)

    # endregion

    # region _filtered_fastq_generator tests
    def test__filtered_fastq_generator(self):
        fw_fastq_string = """@D00611:256:HKTHMBCXX:1:1101:2675:2244 1:N:0:TTAGGC
CCGGTTCATGCCGCCCATGC
+
IIIIGIIIIIIIIIIIIIGI
@D00611:256:HKTHMBCXX:1:1101:2804:2247 1:N:0:TTAGGC
CCACCATTGGTGTGCTGCA
+
GGGGGGGGIGIIIGIIGII
@D00611:256:HKTHMBCXX:1:1101:3177:2193 1:N:0:TTAGGC
GGTCTGTTGCCTGGTCCATC
+
GGIIIIIIIIGIIIIGIIII
@D00611:256:HKTHMBCXX:1:1101:3125:2247 1:N:0:TTAGGC
TCGGCGTCCCCGGCCAGCCA
+
IIIIIIIIIIIIGIIIIIII
@D00611:256:HKTHMBCXX:1:1101:3436:2196 1:N:0:TTAGGC
GGGCCACTAGGGACAGGATGTTTTAGAGCTAGAAATAGCAAGTTA
+
IIGIIIIIGIIIGGIIIGGIGIIIIIIIGIIIIIIIIIGGIGGII
@D00611:256:HKTHMBCXX:1:1101:3709:2197 1:N:0:TTAGGC
CGAATACACTCTGCAGAGCTT
+
IIIIIIIIIIIIIIIIIIIII
"""
        rv_fastq_string = """@D00611:256:HKTHMBCXX:1:1101:2675:2244 2:N:0:TTAGGC
GGGCTCCGGAGGCCATGCC
+
IIIGIIIIIGGGIIIGIII
@D00611:256:HKTHMBCXX:1:1101:2804:2247 2:N:0:TTAGGC
TAATGTTGCTCGAGTTGGA
+
GGGGGGGGGGGIGGGIIGI
@D00611:256:HKTHMBCXX:1:1101:3177:2193 2:N:0:TTAGGC
CATGCTGGATAAGGCGGC
+
IGIIIIIIGGGGIIGGGG
@D00611:256:HKTHMBCXX:1:1101:3125:2247 2:N:0:TTAGGC
TTAACACTCGACGTTTGA
+
IIIIIIGGII<AGGGIIG
@D00611:256:HKTHMBCXX:1:1101:3436:2196 2:N:0:TTAGGC
AAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACG
+
GGIIIIIIIIIIGIIIIGIIGGIGGIIGIGGIGIIGGGGGGI
@D00611:256:HKTHMBCXX:1:1101:3709:2197 2:N:0:TTAGGC
TGCCGGCTGTCGGTATGTCG
+
IIGIIIIIGGIGGGIIGIGG
"""
        fw_fastq_handler = FastqHandler(fw_fastq_string, True)
        rv_fastq_handler = FastqHandler(rv_fastq_string, True)
        counters = {"num_pairs": 0, "num_pairs_passing": 0}

        num_expected_results = 3
        test_num = 0
        filtered_generator = _filtered_fastq_generator(fw_fastq_handler, rv_fastq_handler, 19, 21, 19, counters)
        for curr_fw_record, curr_rv_record in filtered_generator:
            if test_num == 0:
                self.assertEqual("CGGTTCATGCCGCCCATGC", curr_fw_record.sequence)
                self.assertEqual("IIIGIIIIIIIIIIIIIGI", curr_fw_record.quality)
                self.assertEqual("GGGCTCCGGAGGCCATGCC", curr_rv_record.sequence)
                self.assertEqual("IIIGIIIIIGGGIIIGIII", curr_rv_record.quality)
                self.assertEqual(1, counters["num_pairs"])
                self.assertEqual(1, counters["num_pairs_passing"])
            elif test_num == 1:
                self.assertEqual("CCACCATTGGTGTGCTGCA", curr_fw_record.sequence)
                self.assertEqual("GGGGGGGGIGIIIGIIGII", curr_fw_record.quality)
                self.assertEqual("TAATGTTGCTCGAGTTGGA", curr_rv_record.sequence)
                self.assertEqual("GGGGGGGGGGGIGGGIIGI", curr_rv_record.quality)
                self.assertEqual(2, counters["num_pairs"])
                self.assertEqual(2, counters["num_pairs_passing"])
            elif test_num == 2:
                self.assertEqual("AATACACTCTGCAGAGCTT", curr_fw_record.sequence)
                self.assertEqual("IIIIIIIIIIIIIIIIIII", curr_fw_record.quality)
                self.assertEqual("TGCCGGCTGTCGGTATGTC", curr_rv_record.sequence)
                self.assertEqual("IIGIIIIIGGIGGGIIGIG", curr_rv_record.quality)
                self.assertEqual(6, counters["num_pairs"])
                self.assertEqual(3, counters["num_pairs_passing"])
            else:
                self.fail("generator reported too many records: {0}".format(test_num+1))

            test_num += 1

        self.assertEqual(num_expected_results, test_num)
        self.assertEqual(6, counters["num_pairs"])
        self.assertEqual(3, counters["num_pairs_passing"])

    # endregion

    # region _summarize_counts tests
    def test__summarize_counts(self):
        counters = {"num_pairs": 0}
        single_output = _summarize_counts(counters)
        self.assertEqual("num_pairs:0", single_output)

        counters = {"num_pairs": 0, "num_pairs_passing": 1}
        multiple_output = _summarize_counts(counters)
        self.assertEqual("num_pairs:0,num_pairs_passing:1", multiple_output)
    # endregion