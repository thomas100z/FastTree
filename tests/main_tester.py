from classes import *
import unittest


class MainTester(unittest.TestCase):

    def test_main(self):

        parser = aln_parser.AlignmentParser("./resources/fasttree-input.aln")
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

        # # NNI
        for i in range(1): #range(round(nni_round)):
            if joins_amount > 200:
                tp.recompute(tree.active_nodes)

            tree.nearest_neighbor_interchange()
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
