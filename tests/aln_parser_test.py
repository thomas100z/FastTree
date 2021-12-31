import unittest
from classes import aln_parser


class TestAlignmentParser(unittest.TestCase):

    def test_parser(self):
        parser = aln_parser.AlignmentParser("./resources/test-small.aln")
        data = parser.get_data()
        self.assertEqual(len(data), 8)

        algiments = [
            "CGACCGACTGGGTAATTATATTGGCCTTACA",
            "CGTTCCAATGGACGACTGTATTTGTATTATA",
            "CAATGAGAACGATGACCTAACCCAGATTATA",
            "CAATGAGAACGAGGAAGGAAACGAGTCGATA",
            "CGCCCCCCTTACCGAAGGAATTCGTCTGTCA",
            "CGGCCTTCTTAGCGAAGGAACCCAGCTCGCA",
            "CGGCCTTCTTAGCGAAGGAGTTCGCATAACC",
            "CCCTTCAAATGACGAAGGAATCCAACCTTCC"
        ]

        for i, (node, aln) in enumerate(zip(data, algiments)):
            self.assertEqual(node.alignment, aln)
            self.assertEqual(node.name, str(i))


if __name__ == '__main__':
    unittest.main()
