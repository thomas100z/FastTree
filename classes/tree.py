from __future__ import annotations
from queue import PriorityQueue
from .node import Node
from .distances import Distances
import logging
import random
from .total_profile import TotalProfile

logger = logging.getLogger('FastTree')


class Tree:
    root: Node

    def __init__(self, nodes: list[Node], m: int, N: int, L: int, bootstrap=False, bootstrap_round=50) -> None:
        self.nodes = nodes.copy()
        self.active_nodes = nodes.copy()
        self.m = m
        self.N = N
        self.L = L
        self.do_bootstrap = bootstrap
        self.bootstrap_rounds = bootstrap_round
        self.joins = 0

    def to_newick(self) -> str:
        """
        Constructs newick format representation of the tree.
        :return: newick string
        """
        return f"{self.root.newick()};"

    def set_total_profile(self, total_profile: TotalProfile) -> None:
        self.tp = total_profile

    def save(self, path: str) -> None:
        """
        Saves the newick representation of the tree to a file
        :param path: the path of the file to save
        """
        with open(path, 'w') as file:
            file.write(self.to_newick())

    def join_nodes(self, node_1: Node, node_2: Node) -> Node:
        """
        Joines two nodes by creating a parent node and setting the two nodes as its children
        :param node_1: the first node to join
        :param node_2: the second node to join
        :return:
        """
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
                                                        self.active_nodes,
                                                        self.tp) < min_distance and node != joined_node:
                        joined_node.best_known.node = node

        # check the top_hits lists of all other nodes for the removed nodes and replace with active ancestor
        for node in self.active_nodes:
            # Check top-hits
            if node_1 in node.top_hits:
                del node.top_hits[node_1]
            if node_2 in node.top_hits:
                del node.top_hits[node_2]
            if node_1 in node.top_hits or node_2 in node.top_hits and node != joined_node:
                node.top_hits[joined_node] = Distances.neighbor_join_distance(joined_node, node, self.active_nodes,
                                                                              self.tp)

        return joined_node

    def construct_initial_topology(self) -> None:
        """
        Constructs the initial topology of the tree using the top-hits heuristic
        and neighbor-joining criterion
        """
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
                                                                                self.active_nodes, self.tp)
            for node in m_best_known:
                if not node.is_active:
                    m_best_known.remove(node)

            while len(m_best_known) < min(self.m, len(self.active_nodes)):
                temp = best_knows.get()
                if temp.is_active:
                    m_best_known.append(temp)

            # recompute the neighbor joining criterion
            best = None
            least_distance = float('inf')
            for node in m_best_known:
                if node.is_active:
                    distance = Distances.neighbor_join_distance(node, node.best_known.node, self.active_nodes, self.tp)
                    node.best_known.distance = distance
                    if distance < least_distance:
                        best = node
                        least_distance = distance
                else:
                    m_best_known.remove(node)
                    m_best_known.append(best_knows.get())

            # log initial best nodes
            logger.debug(f'{best=}\t, {len(m_best_known)}\t{len(self.active_nodes)}')

            # perform hill climbing for the current two best nodes
            node_1 = best
            node_2 = best.best_known.node

            selected_1 = node_1
            selected_2 = node_2

            # Local Hill Climbing
            for node in node_1.top_hits:
                distance = Distances.neighbor_join_distance(node_1, node, self.active_nodes, self.tp)
                if distance < least_distance:
                    selected_1 = node_1
                    selected_2 = node
                    least_distance = distance

            for node in node_2.top_hits:
                distance = Distances.neighbor_join_distance(node, node_2, self.active_nodes, self.tp)
                if distance < least_distance:
                    selected_1 = node
                    selected_2 = node_2
                    least_distance = distance

            # log final best nodes after hill climbing
            logger.debug(f"joining nodes: {selected_1.name} -  {selected_2.name}")

            joined_node = self.join_nodes(selected_1, selected_2)
            best_knows.put(joined_node)

        # save the last remaining active node as the root of the tree
        self.root = self.active_nodes.pop()
        logger.debug(f'Initial topology:\t{self.to_newick()}')

    def nearest_neighbor_interchange(self) -> None:
        """
        Perform nearest neighbor interchange to evaluate if a split is favorable
        """
        queue = [node for node in self.nodes if node.is_leaf]

        # for all nodes check if they support a split in postorder
        while queue:
            current_node = queue.pop()
            # because all nodes have two children the only condition is to have 3 parents
            if current_node.parent and current_node.parent.parent and current_node.parent.parent.parent:

                # determine the a, b, c and d node for evaluating a different topology
                a = current_node
                b = current_node.get_sibling()
                c = current_node.parent.get_sibling()
                d = a.parent.parent.parent

                logger.debug(f'topology being evaluated: \ta=:{a.name}\tb:{b.name}\t\tc:{c.name}\td:{d.name}')

                # topology abcd
                d_abcd = Distances.log_corrected_profile_distance(a.profile, b.profile) + \
                         Distances.log_corrected_profile_distance(c.profile, d.profile)

                # topology acbd
                d_acbd = Distances.log_corrected_profile_distance(a.profile, c.profile) + \
                         Distances.log_corrected_profile_distance(b.profile, d.profile)

                # topology adbc
                d_bcad = Distances.log_corrected_profile_distance(b.profile, c.profile) + \
                         Distances.log_corrected_profile_distance(a.profile, d.profile)

                logger.debug(f'topology distances - d_abcd:{d_abcd}\t d_acbd:{d_acbd}\t d_bcad:{d_bcad}')

                if d_bcad < min(d_abcd, d_acbd):
                    logger.debug(f'switching nodes: {d.name} - {b.name}')

                    if self.do_bootstrap:
                        if self.bootstrap(a, b, c, d, "d_bcad"):
                            self.switch_nodes(b, a)
                    else:
                        self.switch_nodes(b, a)

                elif d_acbd < min(d_abcd, d_bcad):
                    logger.debug(f'switching nodes: {c.name} - {b.name}')

                    if self.do_bootstrap:
                        if self.bootstrap(a, b, c, d, "d_acbd"):
                            self.switch_nodes(b, c)
                    else:
                        self.switch_nodes(b, c)

            if current_node.parent and current_node.parent not in queue:
                queue.append(current_node.parent)

        # recompute the profile for all internal nodes
        queue = [node.parent for node in self.nodes if node.is_leaf]

        while queue:
            current_node = queue.pop()
            if current_node != self.root and current_node.parent not in queue:
                queue.append(current_node.parent)

    def bootstrap(self, a: Node, b: Node, c: Node, d: Node, split: str) -> bool:
        result = 0

        for i in range(self.bootstrap_rounds):
            random_columns = sorted(random.sample(range(self.L), round(self.L * 0.8)))

            bootstrap_profile_a = a.profile[:, random_columns]
            bootstrap_profile_b = b.profile[:, random_columns]
            bootstrap_profile_c = c.profile[:, random_columns]
            bootstrap_profile_d = d.profile[:, random_columns]

            # topology abcd
            d_abcd = Distances.log_corrected_profile_distance(bootstrap_profile_a, bootstrap_profile_b) + \
                     Distances.log_corrected_profile_distance(bootstrap_profile_c, bootstrap_profile_d)

            # topology acbd
            d_acbd = Distances.log_corrected_profile_distance(bootstrap_profile_a, bootstrap_profile_c) + \
                     Distances.log_corrected_profile_distance(bootstrap_profile_b, bootstrap_profile_d)

            # topology adbc
            d_bcad = Distances.log_corrected_profile_distance(bootstrap_profile_b, bootstrap_profile_c) + \
                     Distances.log_corrected_profile_distance(bootstrap_profile_a, bootstrap_profile_d)

            if split == "d_acbd":
                result += d_acbd < min(d_abcd, d_bcad)

            if split == "d_bcad":
                result += d_bcad < min(d_abcd, d_acbd)

        logger.debug(f'bootstrap support: {result > (self.bootstrap_rounds / 2)}')
        return result > (self.bootstrap_rounds / 2)

    def switch_nodes(self, node_1: Node, node_2: Node) -> None:
        """
        Switches two nodes in the tree. Can be leaves or subtrees
        :param node_1: the first node to switch
        :param node_2: the second node to switch
        """
        parent_1 = node_1.parent
        parent_2 = node_2.parent

        parent_1.children.remove(node_1)
        parent_2.children.remove(node_2)

        parent_2.children.append(node_1)
        node_1.parent = parent_2

        parent_1.children.append(node_2)
        node_2.parent = parent_1

        parent_2.recompute_profile()
        parent_1.recompute_profile()

        # rename the parent nodes
        parent_1.name = parent_1.rename()
        parent_2.name = parent_2.rename()

        logger.debug(f'new topology:\t{self.to_newick()}')

    def set_top_hits(self) -> None:
        """
        Sets a list of top hits for all the active nodes.
        """
        for current_node in self.nodes:
            if not current_node.top_hits:

                node_distances = {}

                for node in self.nodes:
                    if node != current_node:
                        node_distances[node] = Distances.neighbor_join_distance(current_node, node, self.active_nodes,
                                                                                self.tp)

                        # check for best known
                        if node_distances[node] < node.best_known.distance:
                            node.best_known.distance = node_distances[node]
                            node.best_known.node = current_node

                        # check for best known
                        if node_distances[node] < current_node.best_known.distance:
                            current_node.best_known.distance = node_distances[node]
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
                                                                                         self.active_nodes, self.tp)
                                # check for best known
                                if node_distances[node] < node.best_known.distance:
                                    node.best_known.distance = node_distances[node]
                                    node.best_known.node = current_node

                                # check for best known
                                if node_distances[node] < current_node.best_known.distance:
                                    current_node.best_known.distance = node_distances[node]
                                    current_node.best_known.node = node

                        # set m most similar
                        current_node.top_hits = {k: v for i, (k, v) in
                                                 enumerate(sorted(other_distances.items(), key=lambda x: x[1]))
                                                 if i < self.m}

    def set_top_hits_node(self, new_node: Node, children_top_hits: list[Node], children: list[Node]) -> None:
        """
        Sets top hits list for a specific node, that has been created during joining
        :param children: combined children of the joined nodes
        :param new_node: the newly joined node
        :param children_top_hits: list of top hits
        """
        if all(elem not in self.active_nodes for elem in children_top_hits):
            children_top_hits = self.active_nodes

        node_distances: dict[Node, float] = {}
        for node in children_top_hits:
            if node != new_node:
                node_distances[node] = Distances.neighbor_join_distance(new_node, node, self.active_nodes, self.tp)

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
        """
        Calculates the branch length for all nodes
        """
        r = self.root

        for node in self.nodes:

            if node == self.root:
                node.branch_length = 0

            # if current node is leaf node A
            elif node.is_leaf:

                a = node

                b = node.get_sibling()
                node.branch_length = (Distances.log_corrected_profile_distance(a.profile, r.profile) +
                                   Distances.log_corrected_profile_distance(a.profile, b.profile) -
                                   Distances.log_corrected_profile_distance(b.profile, r.profile)) / 2

                if node.branch_length < 0:
                    b.branch_length -= node.branch_length
                    node.branch_length = 0

            # if current node is internal node
            else:
                # current node is AB, consisting of nodes A and B (calculate d(AB,r))
                ab = node
                a = node.children[0]
                b = node.children[1]
                c = node.get_sibling()
                node.branch_length = (Distances.log_corrected_profile_distance(a.profile, r.profile) +
                                    Distances.log_corrected_profile_distance(a.profile, c.profile) +
                                    Distances.log_corrected_profile_distance(b.profile, r.profile) +
                                    Distances.log_corrected_profile_distance(b.profile, c.profile)) / 4 - \
                                   (Distances.log_corrected_profile_distance(a.profile, b.profile) +
                                    Distances.log_corrected_profile_distance(r.profile, c.profile)) / 2

                if node.branch_length < 0:
                    c.branch_length -= node.branch_length
                    node.branch_length = 0

            logger.debug(f'node:{node.name}\tbranch length:{node.branch_length}')
