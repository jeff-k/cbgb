from scipy.sparse.csgraph import laplacian
from numpy.linalg import det
from cbgb import kmerise
import numpy as np

from scipy.special import gammaln

def perplexity4(seq, k=12):
    """in theory the base4-perplexity should suggest the best kmer size(?)
    """
    kmers = {}
    n = 0
    for kmer in kmerise(seq, k=k):
        if kmer not in kmers:
            kmers[kmer] = 0
        kmers[kmer] += 1
        n += 1

    # pow(4, -\sum_{kmer}^{kmers} p_{kmer} log_4 p_{kmer})
    b = np.log(4)
    dist = np.array(list(map(lambda x: (np.log((x/n))/b) * (x/n), kmers.values())))
    return -1 * np.sum(dist)
#    return np.power(4, -1 * np.sum(dist))

def ec(adj):
    """count eulerian cycles w/ BEST theorem
    """
    def degree(v):
        """sum of row vector from adjacency matrix is the out degree"""
        return sum(adj[v]) - adj[v][v]

    # arborescence; determinant of truncated laplacian(?)
    arbors = det(laplacian(adj)[1:,1:])
    # \prod_{v \in G}(d^{+}(v)-1)!
    prod = np.prod([np.math.factorial(degree(v) - 1) for v in range(0,
        len(adj))])

    return arbors * prod

def lnec(adj):
    """ln of count of eulerian cycles, for numerical stability
    """
    def degree(v):
        return sum(adj[v]) - adj[v][v]

    lnarbors = np.log(det(laplacian(adj)[1:,1:]))
    sumln = np.sum([gammaln(degree(v)) for v in range(0, len(adj))])
    return lnarbors + sumln
