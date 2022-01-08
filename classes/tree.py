from __future__ import annotations
from queue import PriorityQueue
from .node import Node
from .distances import Distances


class Tree:
    root: Node

    def __init__(self, nodes: list[Node], m: int, N: int) -> None:
        self.nodes = nodes.copy()
        self.active_nodes = nodes.copy()
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
        childeren =list(set(list(node_1.top_hits.keys()) + list(node_2.top_hits.keys())))

        self.set_top_hits_node(joined_node, childeren)
        self.joins += 1
        self.active_nodes.append(joined_node)
        self.nodes.append(joined_node)

        # remove from active
        self.active_nodes.remove(node_1)
        self.active_nodes.remove(node_2)
        node_1.active = False
        node_2.active = False

        return joined_node

    def construct_initial_topology(self) -> None:

        while len(self.active_nodes) > 1:

            # construct priority queue to find the m best best-know
            best_knows = PriorityQueue()
            any(best_knows.put(node) for node in self.active_nodes)
            m_best_known = [best_knows.get() for _ in range(self.m)]

            # recompute the neighbor joining criterion ...? TODO: does this do any thing? maybe recompute only among these M? that would make a differnce i think
            best = None
            least_distance = float('inf')
            for node in m_best_known:
                distance = Distances.neighbor_join_distance(node, node.best_known.node, self.active_nodes)
                node.best_known.distance = distance
                if distance < least_distance:
                    best = node
                    least_distance = distance

            # perform hill climbing for the best TODO: set best with new node if better join found
            node_1 = best
            node_2 = best.best_known.node

            for node in node_1.top_hits:
                pass

            for node in node_2.top_hits:
                pass

            self.join_nodes(node_1, node_2)

        self.root = self.active_nodes.pop()


    def nearest_neighbor_interchange(self) -> None:
        raise Exception("to be implemented")

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

                current_node.top_hits = {k: v for i, (k, v) in
                                         enumerate(sorted(node_distances.items(), key=lambda x: x[1]))
                                         if i < self.m}

                # compute for other nodes
                for n in current_node.top_hits:
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
                        current_node.top_hits = {k: v for i, (k, v) in
                                                 enumerate(sorted(other_distances.items(), key=lambda x: x[1]))
                                                 if i < self.m}

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

        new_node.top_hits = {k:v for i, (k, v) in
                             enumerate(sorted(node_distances.items(), key=lambda x: x[1]))
                             if i < self.m}

    def calculate_branch_length(self):
        for node in self.nodes:
            if node.is_leaf:
                raise 'to be implemented'
            else:
                raise 'to be implemented'
