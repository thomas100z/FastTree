from __future__ import annotations
from .node import Node
from .distances import Distances


class Tree:
    root: Node

    def __init__(self, nodes: list[Node], m: int, N: int) -> None:
        self.nodes = nodes
        self.active_nodes = nodes
        self.m = m
        self.N = N
        self.joins = 0

    def to_newick(self) -> str:
        return f"({self.root.print_newick()});"

    def deactivate_node(self, node: Node) -> None:
        self.active_nodes.remove(node)

    def save(self, path: str) -> None:
        with open(path, 'w') as file:
            file.write(self.to_newick())

    def join_nodes(self, node_1: Node, node_2: Node) -> Node:

        joined_node = Node(node_1.name + node_2.name, "", Node.join_profiles(node_1.profile, node_2.profile), False)
        joined_node.add_child(node_1)
        joined_node.add_child(node_2)

        self.set_top_hits_node(joined_node, node_1.top_hits + node_2.top_hits)
        self.joins += 1
        self.active_nodes.append(joined_node)
        self.nodes.append(joined_node)

        return joined_node

    def construct_initial_topology(self) -> None:
        raise "to be implemented"

    def nnis(self) -> None:
        raise "to be implemented"

    def set_top_hits(self) -> None:
        """

        """
        for current_node in self.nodes:
            if not current_node.top_hits:

                node_distances = {}

                for node in self.nodes:
                    if node != current_node:
                        node_distances[node] = Distances.neighbor_join_distance(current_node, node, self.active_nodes)

                        # check for best known
                        if node_distances[node] < node.best_known.distance:
                            node.best_known.distance = node_distances[node]
                            node.best_known.node = current_node

                        # check for best known
                        if node_distances[node] < current_node.best_known.distance:
                            node.best_known.distance = node_distances[node]
                            current_node.best_known.node = node

                # take 2m most similar
                node_distances = {k: v for i, (k, v) in
                                  enumerate(sorted(node_distances.items(), key=lambda x: x[1])) if i < 2 * self.m}
                current_node.top_hits = [k for i, (k, v) in
                                         enumerate(sorted(node_distances.items(), key=lambda x: x[1])) if i < self.m]

                # compute for other nodes
                for n in node_distances:
                    current_node = n
                    if not current_node.top_hits:
                        other_distances = {}

                        for node in node_distances:
                            if node != current_node:
                                other_distances[node] = Distances.neighbor_join_distance(current_node, node,
                                                                                         self.active_nodes)
                                # check for best known
                                if node_distances[node] < node.best_known.distance:
                                    node.best_known.distance = node_distances[node]
                                    node.best_known.node = current_node

                                # check for best known
                                if node_distances[node] < current_node.best_known.distance:
                                    node.best_known.distance = node_distances[node]
                                    current_node.best_known.node = node

                        # set m most similar
                        current_node.top_hits = [k for i, (k, v) in
                                                 enumerate(sorted(other_distances.items(), key=lambda x: x[1]))
                                                 if i < self.m]

    def set_top_hits_node(self, new_node: Node, children_top_hits: list[Node]) -> None:
        """

        :param new_node:
        :param children_top_hits:
        """
        node_distances: dict[Node, float] = {}
        for node in children_top_hits:
            if node != new_node:
                node_distances[node] = Distances.neighbor_join_distance(new_node, node, self.active_nodes)

                # check for best known
                if node_distances[node] < node.best_known.distance:
                    node.best_known.distance = node_distances[node]
                    node.best_known.node = new_node

                # check for best known
                if node_distances[node] < new_node.best_known.distance:
                    node.best_known.distance = node_distances[node]
                    new_node.best_known.node = node

        new_node.top_hits = [k for i, (k, v) in
                             enumerate(sorted(node_distances.items(), key=lambda x: x[1]))
                             if i < self.m]

    def calculate_branch_length(self):
        for node in self.nodes:
            if node.is_leaf:
                raise 'to be implemented'
            else:
                raise 'to be implemented'
