# cbgb
Coloured de Bruijn Graph Builder (& other methods for updating graph-edges)

![CBGB OMGFUG](docs/cbgb.jpg)

```python3
from cbgb.omfug import ec
from cbgb import CdB, kmerize, Edge

seq = "pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy"
for k in range(2, 12):
    g = CdB()
    for kmer in kmerize(seq, k=k):
        g.add(kmer, Edge())
    g.circularize()
    graph, _ = g.to_adj()
    print("{}-mer:\t{} walks".format(k, ec(graph)))
    walk = list(g.walk(start=seq[:k-1]))
    start = walk[-1][:-1]
    print('\t', start + ''.join([s[-1] for s in walk[:-1]]))

```

```
2-mer:	1.2772574395192978e+24 walks
	 pointypoinointy-noint anointyes \ ey my anty-po \ oird \ binty 
3-mer:	464380231679999.6 walks
	 pointyointy-noint my eyes \ anointy \ anointy-pointy bird \ o p
4-mer:	71663615999.99995 walks
	 pointyinty-noint my eyes \ anointy \ anointy-pointy bird \ o po
5-mer:	33177599.999999985 walks
	 pointynty-noint my eyes \ anointy \ anointy-pointy bird \ o poi
6-mer:	36863.999999999935 walks
	 pointyty-nointy \ anoint my eyes \ anointy-pointy bird \ o poin
7-mer:	95.99999999999997 walks
	 pointy \ anoint my eyes \ anointyy-nointy-pointy bird \ o point
8-mer:	3.999999999999999 walks
	 pointy \ anoint my eyes \ anointy-nointy bird \ o pointy-pointy
9-mer:	2.0 walks
	 pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy
10-mer:	1.0 walks
	 pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy
11-mer:	1.0 walks
	 pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy
```

### Monoidal edges
Integers for a multigraph, sets for labels

### Eulerian walk
Linear time assembly for well behaved, unlabelled de Bruijn graphs

### What's subdawg?
Extract a directed acyclic word graph
