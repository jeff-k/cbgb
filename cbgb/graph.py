"""sliding slice string of length k
"""
def kmerize(string, k=12):
    kmers = {}
    for i in range(0, len(string) - (k - 1)):
        kmer = string[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = 0
        kmers[kmer] += 1
    return kmers

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
    def __init__(self, label=None, count=1):
        self.count = count
        if label:
            self.labels = set([label,])
        else:
            self.labels = set([])

    def update(self, other):
        self.count += other.count
        self.labels = self.labels.union(other.labels)

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

        self.edges[(left, right)] = edge

    def walk(self):
        path = [] # the path of edges
        return path

    def compress(self):
        pass

    def to_gfa(self):
        pass

    def subdag(self, start, end):
        pass
