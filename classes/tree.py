from __future__ import annotations
from queue import PriorityQueue
from .node import Node
from .distances import Distances
import logging

logger = logging.getLogger('FastTree')

DEBUG = False


class Tree:
    root: Node

    def __init__(self, nodes: list[Node], m: int, N: int) -> None:
        self.nodes = nodes.copy()
        self.active_nodes = nodes.copy()
        self.m = m
        self.N = N
        self.joins = 0

    def to_newick(self) -> str:
        return f"{self.root.print_newick()};"

    def save(self, path: str) -> None:
        with open(path, 'w') as file:
            file.write(self.to_newick())

    def join_nodes(self, node_1: Node, node_2: Node) -> Node:

        joined_node = Node(node_1.name + node_2.name, "", Node.join_profiles(node_1.profile, node_2.profile), False)
        joined_node.add_child(node_1)
        joined_node.add_child(node_2)
        node_1.parent = joined_node
        node_2.parent = joined_node
        node_1.is_active = False
        node_2.is_active = False

        # remove from active
        self.active_nodes.remove(node_1)
        self.active_nodes.remove(node_2)

        children = list(set(list(node_1.top_hits.keys()) + list(node_2.top_hits.keys())))
        self.set_top_hits_node(joined_node, children, joined_node.children)
        self.joins += 1
        self.active_nodes.append(joined_node)
        self.nodes.append(joined_node)

        if len(self.active_nodes) > 1:
            if joined_node.best_known.node not in self.active_nodes:
                min_distance = float('inf')
                for node in self.active_nodes:
                    if Distances.neighbor_join_distance(joined_node, node,
                                                        self.active_nodes) < min_distance and node != joined_node:
                        joined_node.best_known.node = node

        # check the top_hits lists of all other nodes for the removed nodes and replace with active ancestor
        for node in self.active_nodes:
            # Check top-hits
            if node_1 in node.top_hits:
                del node.top_hits[node_1]
            if node_2 in node.top_hits:
                del node.top_hits[node_2]
            if node_1 in node.top_hits or node_2 in node.top_hits and node != joined_node:
                node.top_hits[joined_node] = Distances.neighbor_join_distance(joined_node, node, self.active_nodes)

        return joined_node

    def construct_initial_topology(self) -> None:

        # construct priority queue to find the m best best-know
        best_knows = PriorityQueue()
        any(best_knows.put(node) for node in self.active_nodes)
        m_best_known = [best_knows.get() for _ in range(self.m)]

        while len(self.active_nodes) > 1:

            # Check that all the best joins show to an active node after a join
            for node in self.active_nodes:
                if node.best_known.node not in self.active_nodes:
                    node.best_known.node = node.best_known.node.parent
                    node.best_known.distance = Distances.neighbor_join_distance(node.best_known.node, node,
                                                                                self.active_nodes)
            for node in m_best_known:
                if not node.is_active:
                    m_best_known.remove(node)

            while len(m_best_known) < min(self.m, len(self.active_nodes)):
                temp = best_knows.get()
                if temp.is_active:
                    m_best_known.append(temp)

            """ TO DO: Periodically refresh the top hit lists"""

            # recompute the neighbor joining criterion
            best = None
            least_distance = float('inf')
            for node in m_best_known:
                if node.is_active:
                    distance = Distances.neighbor_join_distance(node, node.best_known.node, self.active_nodes)
                    node.best_known.distance = distance
                    if distance < least_distance:
                        best = node
                        least_distance = distance
                else:
                    m_best_known.remove(node)
                    m_best_known.append(best_knows.get())

            # perform hill climbing for the best
            logger.debug(f'{best=}\t, {len(m_best_known)}\t{len(self.active_nodes)}')

            node_1 = best
            node_2 = best.best_known.node

            selected_1 = node_1
            selected_2 = node_2

            # Local Hill Climbing
            for node in node_1.top_hits:
                distance = Distances.neighbor_join_distance(node_1, node, self.active_nodes)
                if distance < least_distance:
                    selected_1 = node_1
                    selected_2 = node
                    least_distance = distance

            for node in node_2.top_hits:
                distance = Distances.neighbor_join_distance(node, node_2, self.active_nodes)
                if distance < least_distance:
                    selected_1 = node
                    selected_2 = node_2
                    least_distance = distance

            logger.debug(f"joining nodes: {selected_1.name} -  {selected_2.name}")

            joined_node = self.join_nodes(selected_1, selected_2)
            best_knows.put(joined_node)

        self.root = self.active_nodes.pop()

    def nearest_neighbor_interchange(self) -> None:

        # for all nodes check if they support a split
        for current_node in self.nodes:

            # because all nodes have two children the only condition is to have 3 parents
            if current_node.parent and current_node.parent.parent and current_node.parent.parent.parent:

                # determine the a, b, c and d node for evaluating a different topology
                # todo: fill correct nodes here and make sure to do recompute in the right order
                a = None
                b = None
                c = None
                d = None

                logger.debug(f'topology being evaluated: \ta=:{a.name}\tb:{b.name}\tc:{c.name}\td:{d.name}')

                # topology abcd
                d_abcd = Distances.profile_distance(a.profile, d.profile) + Distances.profile_distance(c.profile,
                                                                                                       d.profile)
                # topology acbd
                d_acbd = Distances.profile_distance(a.profile, c.profile) + Distances.profile_distance(b.profile,
                                                                                                       d.profile)
                # topology adbc
                d_adbc = Distances.profile_distance(a.profile, d.profile) + Distances.profile_distance(b.profile,
                                                                                                       c.profile)
                if d_adbc < min(d_abcd, d_acbd):
                    logger.debug(f'switching nodes: {d.name} - {b.name}')
                    self.switch_nodes(d, b)

                elif d_acbd < min(d_abcd, d_adbc):
                    logger.debug(f'switching nodes: {c.name} - {b.name}')
                    self.switch_nodes(b, c)

    @staticmethod
    def switch_nodes(node_1: Node, node_2: Node) -> None:
        parent_1 = node_1.parent
        parent_2 = node_2.parent

        parent_1.children.remove(node_1)
        parent_2.children.remove(node_2)

        parent_2.children.append(node_1)
        node_1.parent = parent_2

        parent_1.children.append(node_2)
        node_2.parent = parent_1

        parent_1.recompute_profile()
        parent_2.recompute_profile()

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

    def set_top_hits_node(self, new_node: Node, children_top_hits: list[Node], children: list[Node]) -> None:
        """

        :param children:
        :param new_node:
        :param children_top_hits:
        """
        if all(elem not in self.active_nodes for elem in children_top_hits):
            children_top_hits = self.active_nodes

        node_distances: dict[Node, float] = {}
        for node in children_top_hits:
            if node != new_node:
                node_distances[node] = Distances.neighbor_join_distance(new_node, node, self.active_nodes)

                # check for best known
                if node_distances[node] < node.best_known.distance:
                    node.best_known.distance = node_distances[node]
                    node.best_known.node = new_node

                # check for best known
                if node_distances[node] < new_node.best_known.distance and node not in children:
                    new_node.best_known.distance = node_distances[node]
                    new_node.best_known.node = node

        new_node.top_hits = {k: v for i, (k, v) in
                             enumerate(sorted(node_distances.items(), key=lambda x: x[1]))
                             if i < self.m}

    def calculate_branch_length(self):
        for node in self.nodes:
            if node.is_leaf:
                raise 'to be implemented'
            else:
                raise 'to be implemented'
