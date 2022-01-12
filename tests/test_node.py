from unittest import TestCase
from classes import *


class TestNode(TestCase):

    def test_pq_for_nodes(self) -> None:
        """
        Test if a priority queue returns the nodes in the correct order
        :return: None
        """
        node_1 = Node("A", "ATCGCG")
        node_2 = Node("B", "ATCGAA")
        node_3 = Node("C", "ATCGGG")
        node_1.best_known.distance = 1
        node_2.best_known.distance = 5
        node_3.best_known.distance = 3
        nodes = [node_1, node_2, node_3]
        pq = PriorityQueue()
        any(pq.put(node) for node in nodes)
        result = []
        while not pq.empty():
            result.append(pq.get().name)

        self.assertEqual(result, ["A", "C", "B"])

    def test_get_sibling(self) -> None:
        A = Node("A", "ATCGCG")
        B = Node("B", "ATCGAA")
        C = Node("C", "ATCGGG")

        A.children.append(B)
        A.children.append(C)
        B.parent = A
        C.parent = A

        self.assertEqual(B.get_sibling().name, 'C')

    def test_rename(self) -> None:
        A = Node("A", "ATCGCG")
        B = Node("B", "ATCGAA")
        C = Node("C", "ATCGGG")
        A.top_hits = {}
        B.top_hits = {}
        C.top_hits = {}

        tree = Tree([A, B, C], 1, 1)

        AB = tree.join_nodes(A, B)
        ABC = tree.join_nodes(AB, C)

        ABC.rename()
        self.assertEqual(ABC.name, "ABC")

    def test_node_comparison(self) -> None:
        A = Node("A", "ATCGCG")
        B = Node("B", "ATCGAA")
        A.best_known.distance = 1
        B.best_known.distance = 5

        self.assertNotEqual(A, B)
        self.assertEqual(A, A)
        self.assertTrue(A < B)
