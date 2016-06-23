# standard libraries
import unittest

# ccbb libraries
from ccbbucsd.utilities.basic_fastq import FastqHandler

# project-specific libraries
from ccbbucsd.malicrispr.construct_counter import _match_and_count_constructs, get_counter_from_names
from ccbbucsd.malicrispr.grna_position_matcher import GrnaPositionMatcher

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region _trim_seq tests
    def test__match_and_count_constructs(self):
        # 1: three; 2: two; 3: two; 4: none; 5: two; 6: one
        fw_fastqs = """@D00611:278:HK55CBCXX:1:1101:1138:2170 1:N:0:ATCACG
ATGCAAGCTCATTGTGAAC
+
GGGGGGGIGGGGGIGAAGG
@D00611:278:HK55CBCXX:1:1101:1375:2236 1:N:0:ATCACG
TTCGGTACGAAACCAGCAC
+
GIIGIIIGGIIIIGIGIGG
@D00611:278:HK55CBCXX:1:1101:1360:2249 1:N:0:ATCACG
TTCGGTACGAAACCAGCAC
+
<<GGGGG<GGGGGGGGIGG
@D00611:278:HK55CBCXX:1:1101:1696:2177 1:N:0:ATCACC
CAACGGGCGTGCCCGGAAA
+
.G.<........<.<...<
@D00611:278:HK55CBCXX:1:1101:1532:2202 1:N:0:ATCACG
TTCGGTACGAAACCAGCAC
+
GIIGGGGGIIGIIIIIIIG
@D00611:278:HK55CBCXX:1:1101:1793:2194 1:N:0:ATCACG
TGTCTGGCCGCGAAGCAGT
+
IIIIIIIIIIIIIIIGIII
"""

        # 1: two rc; 2: none; 3: one rc; 4: none; 5: one rc; 6: two rc
        rv_fastqs = """@D00611:278:HK55CBCXX:1:1101:1138:2170 2:N:0:ATCACG
GTGCGGGTTTCGTACCGAA
+
GGIIGIIGGGGGAGAGGII
@D00611:278:HK55CBCXX:1:1101:1375:2236 2:N:0:ATCACG
TATGCCGGGACTAGAATGG
+
IIIGGGGIIIGGGGGIGII
@D00611:278:HK55CBCXX:1:1101:1360:2249 2:N:0:ATCACG
ACTGCTTCGCGGCCAGACA
+
IGGGGGAGGIGAGGGIGGA
@D00611:278:HK55CBCXX:1:1101:1696:2177 2:N:0:ATCACC
GACACCCAGGGAGCGCGCC
+
.<...<......<......
@D00611:278:HK55CBCXX:1:1101:1532:2202 2:N:0:ATCACG
ACTGCTTCGCGGCCAGACA
+
GIIIIIIAGIIIIIIIIII
@D00611:278:HK55CBCXX:1:1101:1793:2194 2:N:0:ATCACG
GTGCGGGTTTCGTACCGAA
+
GGIIIGIIGIIIIIIIIII
"""

        fw_fastq_handler = FastqHandler(fw_fastqs, True)
        rv_fastq_handler = FastqHandler(rv_fastqs, True)
        grna_names_and_seqs = [("one", "TGTCTGGCCGCGAAGCAGT"),
                               ("two","TTCGGTACGAAACCCGCAC"),  # note one mismatch to fastq seq
                               ("three", "ATGCAAGCTCATTGTGAAC")]
        construct_names = ["three__one", "two__one", "three__two"]
        construct_counts = get_counter_from_names(construct_names)
        grna_matcher = GrnaPositionMatcher(grna_names_and_seqs, 19)
        output_construct_counts, output_summary_counts = _match_and_count_constructs(grna_matcher, construct_counts,
                                                                                     fw_fastq_handler, rv_fastq_handler)
        self.assertEqual(3, len(output_construct_counts))
        self.assertEqual(0, output_construct_counts["three__one"])
        self.assertEqual(2, output_construct_counts["two__one"])
        self.assertEqual(1, output_construct_counts["three__two"])
        # also finds one instance of one__two, but this one is considered unrecognized as it isn't in the list of
        # recognized construct names

        self.assertEqual(5, len(output_summary_counts))
        self.assertEqual(6, output_summary_counts["num_pairs"])
        self.assertEqual(2, output_summary_counts["num_pairs_unrecognized"])  # two, none and none, none
        self.assertEqual(4, output_summary_counts["num_constructs_found"])
        self.assertEqual(1, output_summary_counts["num_constructs_unrecognized"])  # one__two
        self.assertEqual(3, output_summary_counts["num_constructs_recognized"])

