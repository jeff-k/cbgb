"""generate examples of sequences with underlying structure
"""
import random
import cbgb

refs = ["ZFACAB", "ACFAD", "ZFAXAY", "ZAYAX", "FAYX", "ZFAY", "AYZFAX"]

e = {}
e['A'] = lambda: random.choice(['aaa', 'aac', 'aag', 'aat', 'tat'])
e['B'] = lambda: random.choice(['cac', 'ccg', 'cgg'])
e['C'] = lambda: random.choice(['cgt', 'cga', 'cag'])
e['D'] = lambda: random.choice(['ttt', 'tca', 'taa', 'tga', 'tta'])
e['E'] = lambda: random.choice(['gca', 'gac'])
e['F'] = lambda: random.choice(['tgc', 'ttg', 'tag'])

e['X'] = lambda: random.choice([e['B'](), e['C']()])
e['Y'] = lambda: random.choice([e['D'](), e['C']()])

e['Z'] = lambda: random.choice(['a', 'c', 'g', 't', 'at', 'cg', 'tg', 'gg'])

for r in refs:
    print(r)
    seq = ''.join(['^'] + [e[l]() for l in r] + ['$'])
    print(seq)
    g = cbgb.CdB()
    for kmer in cbgb.kmerise(seq):
        g.add(kmer[:-1], cbgb.Edge(data=kmer[-1]))

    print(g.to_dot())
