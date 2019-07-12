import numpy as np

"""sliding slices for strings of length k
"""
def kmerize(seq, k=12):
    for i in range(0, len(seq) - (k - 1)):
        yield seq[i:i+k]

class Edge:
    """prototypical edge type. monoid.
    """
    def __init__(self, data=None):
        """__init__() should return unit edge
        """
        self.data = data

    def update(self, other):
        """monoid operation
        """
        pass

    def __repr__(self):
        if self.data is not None:
            return "{}".format(self.data)
        return "<Edge>"


class Multiedge(Edge):
    """weighted edge for a multigraph
    """
    def __init__(self):
        self.count = 1

    def update(self, other):
        self.count += other.count

    def __repr__(self):
        return "<Edge: {}>".format(self.count)


class LMG(Edge):
    """carrier type for labelled multigraph
    """
    def __init__(self, label=None, count=0):
        self.count = count
        if label:
            self.labels = set([label,])
        else:
            self.labels = set([])

    def update(self, other):
        self.count += other.count
        self.labels = self.labels.union(other.labels)

    def __repr__(self):
        return "<{} {}>".format(self.labels, self.count)


class CdB:
    """main character
    """
    def __init__(self):
        self.edges = {}
        self.nodes = {}
        self.colours = {}

    def add(self, kmer, edge):
        left, right = kmer[:-1], kmer[1:]
        if left not in self.nodes:
            self.nodes[left] = []
        if right not in self.nodes:
            self.nodes[right] = []

        self.nodes[left].append(right)
        if (left, right) not in self.edges:
            self.edges[(left, right)] = edge
        else:
            self.edges[(left, right)].update(edge)

    def walk(self, start=None):
        """eulerian walk!
        """
        path = [] # the path of edges
        node_path = []
        if start is None:
            start = list(self.nodes.keys())[0]

        def visit(node):
            while len(self.nodes[node]) > 0:
                right = self.nodes[node].pop()
                path.append(self.edges[(node, right)])
                if right not in self.nodes:
                    print("broken loop; not eulerian cycle?")
                    break
                visit(right)
            node_path.append(node)

        visit(start)
        return reversed(node_path)

    def compress(self):
        """concatenate nodes which have outdegree of one with their neighbours
        """
        new = {}
        visited = set()
        for node in self.nodes:
            if node in visited:
                continue
            visited.add(node)
            nn = node
            while len(self.nodes[node]) == 1:
                nn += self.nodes[node][0][-1]
                node = self.nodes[node][0]
                visited.add(node)
            new[nn] = self.nodes[node]

        # proof of concept!! inefficient.
        #new_names = {}
        #for node in new:
        #    for out in new[node]:
        #        for name in new.keys():
        #            if out == name[:len(out)]:
        #                new_names[out] = name

        #for node in new:
        #    new[node] = [new_names[out] for out in new[node]]
        return new

    def to_adj(self):
        """dump graph as adjacency matrix
        """
        node_order = list(sorted(self.nodes.keys()))
        rows = []

        for node in node_order:
            row = [0 for _ in range(0, len(node_order))]
            for out in self.nodes[node]:
                row[node_order.index(out)] += 1
            rows.append(row)

        return np.array(rows), node_order

    def to_gfa(self, gfa_file, csv_file):
        gfa_fd = open(gfa_file, 'w')
        csv_fd = open(csv_file, 'w')
        o = self.compress()
        #o = self.nodes 
        names = list(sorted(o.keys()))
        # every node is an alignment
        print("H\tVN:Z:cbgb-omfug", file=gfa_fd)
        print("Name, Label", file=csv_fd)
        for node in o:
            print("S\t{}\t{}".format(names.index(node), node), file=gfa_fd)
            print("{}, \"{}\"".format(names.index(node), node), file=csv_fd)
            for out in o[node]:
                print("L\t{}\t+\t{}\t+\t1M".format(names.index(node),
                    names.index(out)), file=gfa_fd)


    def circularize(self, edge=None):
        """if the graph was generated from a linear sequence there should be
        exactly two nodes with mismatched in/out degrees. join them up.
        """
        if edge is None:
            edge = Edge()
        ins = {}
        wrong_in = []
        wrong_out = []
        for node in self.nodes:
            if node not in ins:
                ins[node] = []
            for out in self.nodes[node]:
                if out not in ins:
                    ins[out] = []
                ins[out].append(node)

        for node in self.nodes:
            if len(ins[node]) > len(self.nodes[node]):
                wrong_out.append(node)
            if len(ins[node]) < len(self.nodes[node]):
                wrong_in.append(node)

        # two wrongs make a right
        if not (len(wrong_in) == 1 and len(wrong_out) == 1):
            return
        left, right = wrong_in[0], wrong_out[0]
        # what if we add a special start symbol instead

        if right not in self.nodes:
            self.nodes[right] = []
        self.nodes[right].append("$")
        self.nodes["$"] = [left]
        self.edges[(right, "$")] = Edge() # is start the unit edge?
        self.edges[("$", left)] = Edge()

    def subdawg(self, start, end):
        """walk from start to end, breadth first

        idea: remove all in-edges to start and out edges from end, add edges
        between start and end with the correct degree. Start walk from start to
        end. Use it to contruct partial order alignments from
        """

        # remove nodes entering start and leaving end
        for (left, right) in self.edges.keys():
            if right is start or left is end:
                self.edges.pop((left, right), None)
        self.edges[(start, end)] = None

    def path_divergence(self, other):
        pass

    def colourize(self, colour):
        pass

def colourUnion(g1, g2):
    g = CdB()
    for k in g1.nodes:
        g.add(k, LMG(label='red'))
    for k in g2.nodes:
        g.add(k, LMG(label='blue'))

    return g
