from .node import Node
import numpy as np


class TotalProfile:
    active: int
    total_profile: np.array

    def __init__(self, nodes: list[Node]) -> None:
        self.nodes = nodes
        self.recompute()

    def recompute(self) -> None:
        profile = 0
        active = 0

        for node in self.nodes:
            if node.active:
                profile += node.profile
                active += 1

        self.active = active
        self.total_profile = profile / active if active != 0 else 0

