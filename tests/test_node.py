from unittest import TestCase
from classes import *


class TestNode(TestCase):

    def test_pq_for_nodes(self):
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
