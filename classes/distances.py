import math
import numpy as np
from .node import Node
from .total_profile import TotalProfile


class Distances:

    @staticmethod
    def uncorrected_distance(sequence_1: str, sequence_2: str) -> float:
        """
        du: uncorrected sequence distance
        :param sequence_1: the first sequence
        :param sequence_2: the second sequence
        """
        result = 0
        gaps = 0
        for x, y in zip(sequence_1, sequence_2):
            if x != '-' and y != '-':
                result += (x != y)
            else:
                gaps += 1

        return result / (len(sequence_1) - gaps)

    @staticmethod
    def corrected_distance(sequence_1: str, sequence_2: str) -> float:
        """
        d: log-corrected sequence distance
        :param sequence_1: the first sequence
        :param sequence_2: the second sequence
        """
        du = 1 - (4 / 3) * Distances.uncorrected_distance(sequence_1, sequence_2)
        if 0 < du:
            return -3 / 4 * math.log10(du)
        else:
            return 1.5

    @staticmethod
    def profile_distance(profile_1: np.array, profile_2: np.array) -> float:
        """
        D: profile distance
        :param profile_1: the first profile
        :param profile_2: the second profile
        """
        S = 0
        for i in range(4):
            for j in range(len(profile_1[0])):
                S += abs(profile_1[i, j] - profile_2[i, j])

        return S / (4 * len(profile_1[0]))

    @staticmethod
    def up_distance(node: Node) -> float:
        """
        u: up-distance
        u(i) = 0 (Leaf Node)
        u(ij) = D(ij)/2 (Half the profile distance)
        :param node: the given node
        """
        if node.is_leaf:
            return 0.0
        else:
            return Distances.profile_distance(node.children[0].profile, node.children[-1].profile) / 2

    @staticmethod
    def node_distance(node1: Node, node2: Node) -> float:
        return Distances.profile_distance(node1.profile, node2.profile) - Distances.up_distance(
            node1) - Distances.up_distance(node2)

    @staticmethod
    def out_distance(node: Node, active: list) -> float:
        """
        r(i): Out Distance
        :param node: the given node
        :param active: the list of all active nodes
        """
        S = 0
        for n in active:
            if n != node:
                S += Distances.node_distance(node, n)
        return S / (len(active) - 2) if len(active) > 2 else S / (len(active))

    @staticmethod
    def neighbor_join_distance(node1: Node, node2: Node, active: list[Node], total_profile: TotalProfile) -> float:
        """
        Neighbor-Joining Distance: du(i,j) - r(i) - r(j)
        :param total_profile: total profile of all the nodes
        :param node1: the first node
        :param node2: the second node
        :param active: the list of all active nodes
        """
        return Distances.node_distance(node1, node2) - Distances.total_profile_out_distance(node1, total_profile,
                                                                                                    active) - \
                       Distances.total_profile_out_distance(node2, total_profile, active)

    @staticmethod
    def total_profile_out_distance(node: Node, tp: TotalProfile, active: list[Node]) -> float:
        """
        r(i): Out Distance = n∆(i, T) − ∆(i, i) − (n − 1)u(i) + u(i) - Σ u(j)
        :param node: the given node
        :param tp: the total profile of the tree
        :param active: the list of all active nodes
        """
        S = 0
        for n in active:
            S += Distances.up_distance(node)

        result = len(active) * Distances.profile_distance(node.profile, tp.total_profile) - \
                 Distances.average_node_children_distance(node) - (len(active) - 2) * Distances.up_distance(node) - S

        return result / (len(active) - 2) if len(active) > 2 else result / (len(active))

    @staticmethod
    def average_node_children_distance(node: Node) -> float:
        """
        ∆(i, i): the average distance between children of i, including self-comparisons
        :param node: the given node under examination
        """
        if len(node.children) == 0: return 0
        S = 0
        for child in node.children:
            for child2 in node.children:
                S += Distances.profile_distance(child.profile, child2.profile)
        return S / (len(node.children) * len(node.children))

    @staticmethod
    def log_corrected_profile_distance(profile_1: np.array, profile_2: np.array) -> float:
        """
        Calculates the log corrected profile distance between two profiles: - 3/4 log( 1 - 4/3 d)
        where d is the profile distance
        :param profile_1: the first profile
        :param profile_2: the second profile
        :return:
        """
        du = 1 - (4 / 3) * Distances.profile_distance(profile_1, profile_2)
        if 0 < du:
            return -3 / 4 * math.log10(du)
        else:
            return 1.5