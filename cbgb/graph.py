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
        pass

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

    def to_gfa(self):
        pass

    def circularize(self, edge=Edge()):
        """if the graph was generated from a linear sequence there should be
        exactly two nodes with mismatched in/out degrees. join them up.
        """

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
        if right not in self.nodes:
            self.nodes[right] = []
        self.nodes[right].append(left)
        self.edges[(right, left)] = edge

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
