from .node import Node
import numpy as np


class TotalProfile:
    active: int
    total_profile: np.array

    def __init__(self, nodes: list) -> None:
        self.nodes = nodes
        self.recompute(nodes)

    def recompute(self, active_nodes: list[Node]) -> None:
        """
        recomputes the total profile from all the active nodes
        :param active_nodes: list of active nodes
        """
        self.nodes = active_nodes
        profile = 0
        active = 0

        for node in self.nodes:
            profile += node.profile
            active += 1

        self.active = active
        self.total_profile = profile / active if active != 0 else 0

