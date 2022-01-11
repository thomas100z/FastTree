from __future__ import annotations
from functools import total_ordering
import numpy as np
import logging

logger = logging.getLogger('FastTree')


class BestKnown:
    node: Node
    distance: float

    def __init__(self) -> None:
        self.distance = float('inf')


@total_ordering
class Node:

    def __init__(self, name: str, alignment: str, profile: np.array = None, is_leaf: bool = True) -> None:
        self.children = []
        self.parent = None
        self.alignment = alignment
        self.name = name
        self.is_leaf = is_leaf
        self.is_active = True
        self.branch_length = 1
        self.best_known = BestKnown()

        if is_leaf:
            self.profile = self.form_profile()
        else:
            self.profile = profile

        self.top_hits = []

    def _is_valid_operand(self, other):
        return hasattr(other, "best_known") and hasattr(self, "best_known")

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.best_known.distance < other.best_known.distance

    def __hash__(self):
        return hash(self.name)

    def add_child(self, node: Node) -> None:
        self.children.append(node)

    def print_newick(self) -> str:
        result = ""
        if self.is_leaf:
            result += f'{self.name}:{self.branch_length}'
        else:
            self.children.sort(key=lambda x: x.is_leaf, reverse=True)
            result += f'({",".join([child.print_newick() for child in self.children])})'

        return result

    @staticmethod
    def join_profiles(profile1: np.array, profile2: np.array) -> np.array:
        return np.mean(np.array([profile1, profile2]), axis=0)

    def form_profile(self, psuesdocount=False) -> np.array[float]:
        sequences = [self.alignment]
        base_value = 4 if psuesdocount else 0

        bases = ['A', 'C', 'G', 'T']

        if psuesdocount:
            profile = np.ones(shape=(4, len(sequences[0])))
        else:
            profile = np.zeros(shape=(4, len(sequences[0])))

        for i in range(0, len(sequences)):
            for j in range(0, len(sequences[i])):
                if sequences[i][j] == " " or sequences[i][j] == "-":
                    continue
                profile[bases.index(sequences[i][j])][j] += 1

        return profile / (len(sequences) + base_value)

    def recompute_profile(self) -> None:
        self.profile = np.mean(np.array([node.profile for node in self.children]), axis=0)

    def get_sibling(self) -> Node:
        parent_children = self.parent.children
        return parent_children[parent_children.index(self) - 1]

