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


    def test_tree_topology(self) -> None:
        # testing some basic tree structures
        A = Node("A", "ATCGCG")
        B = Node("B", "ATCGAA")
        C = Node("C", "ATCGGG")

        tree = Tree([A,B,C], 2,3)
        tree.set_top_hits()
        tree.construct_initial_topology()

        self.assertEqual(tree.root.name, "ABC")
        self.assertEqual([i.name for i in tree.root.children], ['C', 'AB'])
        self.assertEqual(A.get_sibling(), B)


    def test_switch_node(self) -> None:

        A = Node("A", "ATCGCG")
        B = Node("B", "ATCGAA")
        C = Node("C", "ATCGGG")

        tree = Tree([A,B,C], 2,3)
        tree.set_top_hits()
        tree.construct_initial_topology()

        tree.switch_nodes(A,C)

        self.assertEqual(C.get_sibling(), B)
        self.assertEqual(A.parent, tree.root)
        self.assertEqual(C, B.get_sibling())
        self.assertEqual(C.parent , A.get_sibling())