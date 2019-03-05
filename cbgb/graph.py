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
        self.data = data

    def update(self, other):
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

        self.nodes[left].append(right)

        if (left, right) not in self.edges:
            self.edges[(left, right)] = edge
        else:
            self.edges[(left, right)].update(edge)

    def walk(self, start=None):
        """eulerian walk!
        """
        path = [] # the path of edges

        if start is None:
            start = list(self.nodes.keys())[0]

        def visit(node):
            while len(self.nodes[node]) > 0:
                right = self.nodes[node].pop()
                path.append(self.edges[(node, right)])
                if right not in self.nodes:
                    print("broken loop; not eulerian cycle?")
                    break # if there's a broken loop
                visit(right)
        visit(start)
        return path

    def compress(self):
        pass

    def to_adj(self):
        """dump graph as adjacency matrix
        """
        return np.array([[]])

    def to_gfa(self):
        pass

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

        # create contrived seam
        for kmer in kmerize(start + end):
            self.add(kmer, Edge(data=kmer[-1]))
        pass

def test(seq1, seq2):
    g = CdB()
    for kmer in kmerize(seq1):
        g.add(kmer, LMG("seq1", 1))
    for kmer in kmerize(seq2):
        g.add(kmer, LMG("seq2", 1))

    return g.edges

def assemble(kmers):
    c = CdB()
    for kmer in kmers:
        c.add(kmer, Edge(data=kmer[-1]))
    return list(map(str, c.walk()))

if __name__ == "__main__":
    #g = test("ATCGACTACGACTGCTGCATCGATCGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTACGATCGATCGATCG",
    #        "ATCTACTTTTTTTTTTTTTTCGCGGCAGCGATCATCGATCGACGACGATACGATCGATCGACTCGATCGATCGACCGGGGGGGGGGGGGGGACACTCGATCGATCGATCGATCACTACGATGCATGCATGCAT")
    #for k in g:
    #    print(g[k])
    #print(g)

    c = CdB()
    rando = "TAGTGAATCGAAGCGCGGCTTCAGAATACCGTTTTGGCTACCTGATACAAAGCCCATCGTGGTCCTCAGATATCGTGCACGTAGAGTTGCACCGCACGCATGTGGAATTAGTGGCGAAGTACGATTCCAAGACCGACGTACGATACAACTATGCGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCT"
    for kmer in kmerize(rando):
        c.add(kmer, Edge(data=kmer[-1]))
    c.subdawg(rando[:12], rando[-12:])
    print(''.join(map(str, c.walk())))
    print(''.join(assemble(kmerize("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))))
    randoloop = "TAGTGAATCGAAGCGCGGCTTCAGAATACCGTTTTGGCTACCTGATACAAAGCCCATCGTGGTCCTCAGATATCGTGCACGTAGAGTTGCACCGCACGCATGTGGAATTAGTGGCGAAGTACGATTCCAAGACCGACGTACGATACAACTATGCGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCTTAGTGAATCGAAGCGCGGC"
    print(''.join(assemble(kmerize(randoloop))))
