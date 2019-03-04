"""sliding slices for strings of length k
"""
def kmerize(seq, k=12):
    for i in range(0, len(seq) - (k - 1)):
        yield seq[i:i+k]

class Edge:
    """prototypical edge type. monoid.
    """
    def __init__(self):
        pass

    def update(self, other):
        return False

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
        self.nodes = set([])

    def add(self, kmer, edge):
        left, right = kmer[1:], kmer[:-1]
        self.nodes.add(left)
        self.nodes.add(right)

        if (left, right) not in self.edges:
            self.edges[(left, right)] = LMG()
        self.edges[(left, right)].update(edge)

    def walk(self):
        """eulerian walk!
        """
        path = [] # the path of edges
        return path

    def compress(self):
        pass

    def to_gfa(self):
        pass

    def subdag(self, start, end):
        """walk from start to end, breadth first
        """
        pass

def test(seq1, seq2):
    g = CdB()
    for kmer in kmerize(seq1):
        g.add(kmer, LMG("seq1", 1))
    for kmer in kmerize(seq2):
        g.add(kmer, LMG("seq2", 1))

    return g.edges

if __name__ == "__main__":
    g = test("ATCGACTACGACTGCTGCATCGATCGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTACGATCGATCGATCG",
            "ATCTACTTTTTTTTTTTTTTCGCGGCAGCGATCATCGATCGACGACGATACGATCGATCGACTCGATCGATCGACCGGGGGGGGGGGGGGGACACTCGATCGATCGATCGATCACTACGATGCATGCATGCAT")
    for k in g:
        print(g[k])
    print(g)
