from typing import Generator, Any
import numpy as np

"""sliding slices for strings of length k
"""

Seq = str
Monoid = Any

def kmerise(seq: Seq, k: int = 12) -> Generator[Seq, None, None]:
    for i in range(0, len(seq) - (k - 1)):
        yield seq[i:i+k]

class Edge:
    """example of monoidal edge; must implement addition and unit
    """
    def __init__(self, data: Monoid = None):
        """__init__() should return unit edge
        """
        self.data = data
        self.value = 1

    def __add__(self, other: Edge) -> Edge:
        """monoid operation
        """
        data = None
        if self.data:
            data = self.data
        elif other.data:
            data = self.data

        e = Edge(data=data)
        e.value = self.value + other.value
        return e

    def __repr__(self) -> str:
        if self.data is not None:
            return f"<{self.data}: {self.value}>"
        return str(self.value)


class LMG:
    """labelled multigraph (coloured edges)
    """
    def __init__(self):
        self.counts = {}

    def __add__(self, other: LMG) -> LMG:
        for k in other.counts:
            if k in self.counts:
                self.counts[k] += other.counts[k]
            else:
                self.counts[k] = other.counts[k]
        return self

    def __repr__(self) -> str:
        s = ', '.join(["f{k}: {self.counts[k]}" for k in self.counts])
        return f"<{s}>"


def label(graph: LMG, colour: Any) -> None:
    """colour all edges of a graph the same way
    """
    any(map(lambda e: LMG(e, label=colour), graph.nodes.values))


class CdB:
    """main character
    """
    def __init__(self, kmers: Generator[Seq, None, None] | None):
        self.nodes = {}

        if kmers:
            map(self.add, kmers)

    def add(self, kmer: str, edge: Edge):
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

    def walk(self, start: Seq | None) -> Generator[Seq, None, None]:
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

    def compress(self) -> GenGraph:
        """concatenate nodes which have outdegree of one with their neighbours,
        yielding a genome reference graph
        """
        new = GenGraph()
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

    def to_dot(self) -> str:
        dot = ['digraph G {']
        for v in self.nodes:
            for e in self.nodes[v]:
                n = v[1:] + e
                dot.append(f'  "{v}" -> "{n}" [label="{e}"];')
        dot.append('}')
        return '\n'.join(dot)

    def circularise(self, edge=None):
        """if the graph was generated from a linear sequence there should be
        two riversides in your city with an odd number of bridges. join them up.
        """
        if edge is None:
            edge = 1
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
        end. Use it to contruct partial order alignments
        """

        # remove nodes entering start and leaving end
        for (left, right) in self.edges.keys():
            if right is start or left is end:
                self.edges.pop((left, right), None)
        self.edges[(start, end)] = None

    def path_divergence(self, other):
        pass

class GenGraph:
    """genome reference graph
    """
    def __init__(self):
        self.nodes = {}
        self.labels = {}

    def to_gfa(self, gfa_fd=None):
        names = list(sorted(self.nodes))

        # every node is an alignment
        print("H\tVN:Z:cbgb-omfug", file=gfa_fd)
        for node in self.nodes:
            print(f"S\t{names.index(node)}\t{node}", file=gfa_fd)
            for out in self.nodes[node]:
                print(f"L\t{names.index(node)}\t+\t{names.index(out)}\t+\t1M",
                        file=gfa_fd)

    def gfa_csv(self, csv_fd=None):
        print("Name, Label", file=csv_fd)
        for node in self.nodes:
            print(f'{names.index(node)}, "{node}"', file=csv_fd)

    def to_dot(self) -> str:
        dot = ['digraph G {']
        for v in self.nodes:
            for e in self.nodes[v]:
                n = v[1:] + e
                dot.append(f'  "{v}" -> "{n}" [label="{e}"];')
        dot.append('}')
        return '\n'.join(dot)

    def kmers(self):
        """return kmers of a compressed dBG
        """
        return []
