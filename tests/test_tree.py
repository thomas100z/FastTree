from unittest import TestCase
from classes import *


class TestTree(TestCase):

    def test_join_nodes(self):
        node_1 = Node("A", "ATCGCG")
        node_2 = Node("C", "ATCGAA")
        nodes = [node_1, node_2]
        N = len(nodes)
        m = round(math.sqrt(N))
        tree = Tree(nodes, m, N)
        tree.set_top_hits()
        new_node = tree.join_nodes(node_1, node_2)

        # assert name join
        self.assertEqual(new_node.name, "AC")
        # assert top hits set
        self.assertTrue(len(new_node.top_hits) > 0)


