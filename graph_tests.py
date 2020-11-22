from cbgb.omfug import ec, perplexity4
from cbgb import CdB, kmerize, LMG, Edge

def colour(seq1, seq2):
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
    g = colour("ATCGACTACGACTGCTGCATCGATCGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTACGATCGATCGATCG",
            "ATCTACTTTTTTTTTTTTTTCGCGGCAGCGATCATCGATCGACGACGATACGATCGATCGACTCGATCGATCGACCGGGGGGGGGGGGGGGACACTCGATCGATCGATCGATCACTACGATGCATGCATGCAT")
    #for k in g:
    #    print(g[k])
    #print(g)

    c = CdB()
    rando = "TAGTGAATCGAAGCGCGGCTTCAGAATACCGTTTTGGCTACCTGATACAAAGCCCATCGTGGTCCTCAGATATCGTGCACGTAGAGTTGCACCGCACGCATGTGGAATTAGTGGCGAAGTACGATTCCAAGACCGACGTACGATACAACTATGCGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCT"
    for kmer in kmerize(rando):
        c.add(kmer, Edge(data=kmer[-1]))
    c.circularize()
    c.subdawg(rando[:12], rando[-12:])
    #print(''.join(map(str, c.walk())))
    #print(''.join(assemble(kmerize("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))))
    randoloop = "TAGTGAATCGAAGCGCGGCTTCAGAATACCGTTTTGGCTACCTGATACAAAGCCCATCGTGGTCCTCAGATATCGTGCACGTAGAGTTGCACCGCACGCATGTGGAATTAGTGGCGAAGTACGATTCCAAGACCGACGTACGATACAACTATGCGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCTTAGTGAATCGAAGCGCGGC"
    #print(''.join(assemble(kmerize(randoloop))))

    x = CdB()
    for kmer in kmerize("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", k=3):
        x.add(kmer, Edge(data=kmer[-1]))
    x.circularize()
    graph, labels = x.to_adj()
    print(graph)
    print(ec(graph))
    print(''.join(assemble(kmerize("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTt", k=3))))

    t = CdB()
    test2 = "ZABCDABEFABY"
    print(test2)
    for kmer in kmerize(test2, k=3):
        t.add(kmer, Edge(data=kmer[-1]))
    graph1, labels = t.to_adj()
    print(graph1, labels)
    t.circularize()
    graph2, labels = t.to_adj()
    print(graph2)
    print(ec(graph2))
    print(''.join(assemble(kmerize(test2, k=3))))

    print(x.to_dot())
