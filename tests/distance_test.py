import unittest
from classes import *


class DistanceTest(unittest.TestCase):
    parser = aln_parser.AlignmentParser("./resources/test-small.aln")
    nodes = parser.get_data()
    DEBUG = False

    def test_total_profile(self):
        pass
        # if self.DEBUG:
        #     print('test: asserting O(n^2) distance between all nodes is same as total profile distance:')

        tp = TotalProfile(self.nodes)
        #
        # for node in self.nodes:
        #     tp_distance = Distances.profile_distance(node.profile, tp.total_profile)
        #     kw_distance = sum([Distances.uncorrected_distance(node.alignment, n.alignment) for n in self.nodes if n!= node])
        #     self.assertEqual(tp_distance, kw_distance)

        # for node in self.nodes:
        #     d = 1000
        #     name = ''
        #     for other in self.nodes:
        #         if node != other:
        #             distance = Distances.uncorrected_distance(node.alignment, other.alignment)
        #             if distance < d:
        #                 d =  distance
        #                 name = other.name

        for node in self.nodes:
            print(f'profile distance from node {node.name} to total profile {Distances.profile_distance(node.profile, tp.total_profile)}')

            #print(node.name, d, name )

if __name__ == '__main__':
    unittest.main()
