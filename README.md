# cbgb
Coloured de Bruijn Graph Builder (& other methods for updating graph-edges)

This is a reference implementation for exploring graph algorithms on small genome de Bruin graph assemblies. There is no intention to support concise datastructures.

![CBGB OMGFUG](docs/cbgb.jpg)

```python3
from cbgb.omfug import ec
from cbgb import CdB, kmerize, Edge

seq = "pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy"
for k in range(2, 11):
    g = CdB()
    for kmer in kmerize(seq, k=k):
        g.add(kmer, Edge())
    g.circularize()
    graph, _ = g.to_adj()
    print("{}-mer:\te^{} cycles".format(k, lnec(graph)))
    walk = list(g.walk(start=seq[:k-1]))
    start = walk[-1][:-1]
    print('\t', start + ''.join([s[k-2] for s in walk[:-1]]))
```
Which yields
```
2-mer:	e^62.73851624705093 cycles
	 pointypoinointy-noint anointyes \ ey my anty-po \ oird \ binty 
3-mer:	e^33.7717247974044 cycles
	 pointyointy-noint my eyes \ anointy \ anointy-pointy bird \ o p
4-mer:	e^24.99524900805808 cycles
	 pointyinty-noint my eyes \ anointy \ anointy-pointy bird \ o po
5-mer:	e^17.317385507379875 cycles
	 pointynty-noint my eyes \ anointy \ anointy-pointy bird \ o poi
6-mer:	e^10.514990744055561 cycles
	 pointyty-nointy \ anoint my eyes \ anointy-pointy bird \ o poin
7-mer:	e^4.564348191467836 cycles
	 pointy \ anoint my eyes \ anointyy-nointy-pointy bird \ o point
8-mer:	e^1.3862943611198904 cycles
	 pointy \ anoint my eyes \ anointy-nointy bird \ o pointy-pointy
9-mer:	e^0.6931471805599453 cycles
	 pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy
10-mer:	e^0.0 cycles
	 pointy bird \ o pointy-pointy \ anoint my eyes \ anointy-nointy
```
Using the `compress()` and `to_gfa()` methods, we can visualize the output with as a GFA file:
![Pointy bird alignment for k=9](docs/mergedbird.png)

### Monoidal edges
Integers for a multigraph, sets for labels

### Eulerian walk
Linear time assembly for well behaved, unlabelled de Bruijn graphs

### What's subdawg?
Extract a directed acyclic word subgraph from between two nodes

### Other methods
  * Enumerate Eulerian cycles with BEST theorem (efficiently as ln \left | L_G^* \right | + \sum_{u \in V} ln\Gamma d(u))
  * estimate optimal k-mer sizes with entropy/perplexity
