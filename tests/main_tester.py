from classes import *
import unittest


class MainTester(unittest.TestCase):

    def test_main(self):

        parser = aln_parser.AlignmentParser("./resources/test-small.aln")
        nodes = parser.get_data()


        print('\nAll nodes:')
        for node in nodes:
            print(f'{node.name=}\t{node.alignment}\t {node.profile=}')

        print('total profile:')
        tp = TotalProfile(nodes)
        print(f'{tp.total_profile=}')

        # set top hits list
        N = len(nodes)
        m = round(math.sqrt(N))
        tree = Tree(nodes, m, N)
        a = 4
        joins_amount = 0
        nni_round = math.log(N) / math.log(2) + 1
        print(f'{m=}\t{N=}')

        # set top hist list for every node
        tree.set_top_hits()

        for node in nodes:
            print(f'{node.name=} top hits: {[_.name for _ in node.top_hits]}\t {node.best_known.distance=} {node.best_known.node.name=}')

        # set the best join
        tree.construct_initial_topology()


        for node in nodes:
            print(f'{node=}\t{node.is_active=}')
        print(tree.to_newick())
        # # NNI
        previous_changes = 0
        for i in range(round(nni_round)):
            previous_changes = tree.joins

            if tree.joins > 200:
                tp.recompute(tree.active_nodes)

            tree.nearest_neighbor_interchange()

            if previous_changes == tree.joins:
                break

        #
        # # local bootstrap
        #
        # # branch length
        # for node in tree.nodes:
        #     pass
        # # print the tree
        print(tree.to_newick())


if __name__ == '__main__':
    unittest.main()
