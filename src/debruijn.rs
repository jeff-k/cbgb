use core::hash::Hash;
use std::cmp::Eq;
use std::collections::HashMap;
//use std::collections::HashSet;
use std::fmt::{Debug, Display};

use bio_seq::prelude::*;

use crate::Monoid;

pub struct DeBruijn<A: Codec + Eq, const K: usize, Edge: Hash + Monoid> {
    map: HashMap<Kmer<A, K>, (A, Edge)>,
}

impl<A: Codec + Eq + Display + Debug, const K: usize, Edge: Monoid + Hash + Display>
    DeBruijn<A, K, Edge>
{
    /*
    fn compress(self) -> GenomeGraph<A, Edge> {
        GenomeGraph {
            edges: HashMap::new(),
            vertices: HashSet::new(),
        }
    }
    */

    fn from_seq(seq: &SeqSlice<A>) -> Self {
        if seq.len() < K {
            panic!("Error: seq too short");
        }
        let mut kmers = HashMap::new();
        for subseq in seq.windows(K + 1) {
            let kmer: Kmer<A, K> = subseq[..K].into();
            let out: A = A::unsafe_from_bits(u8::from(&subseq[K]));
            let (_v, ref mut e) = kmers.entry(kmer).or_insert((out, Edge::ZERO));
            *e = e.addm(&Edge::ONE);
        }
        DeBruijn { map: kmers }
    }

    /*
    fn push_kmer(mut self, kmer: Kmer<A, N>) {
        self.map.insert(kmer, (Edge::zero());
    }
    */
    fn dump(self) {
        for (i, (kmer, (out, edge))) in self.map.into_iter().enumerate() {
            println!("{i}\t{kmer}\t\t{out}\t{edge}");
        }
    }
}

/*
impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphBase for DeBruijn<A, N, Edge> {
    type EdgeId = Edge;
    type NodeId = Kmer<A, N>;
}
*/

/*
impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphRef for &DeBruijn<A, N, Edge>
where [(); N+1]:,{}
*/

/*
pub struct NeighbourIter<A: Codec, const N: usize> {
    mers: HashMap<Kmer<A, N>, usize>,
}

impl<A: Codec, const N: usize> Iterator for NeighbourIter<A, N> {
    type Item = Kmer<A, N>;
    fn next(&mut self) -> Option<Kmer<A, N>> {
        None
    }
}

impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> IntoNeighbors
    for &DeBruijn<A, N, Edge>
{
    type Neighbors = NeighbourIter<A, N>;

    fn neighbors(self, node: <Self as GraphBase>::NodeId) -> Self::Neighbors {
        todo!()
    }
}
*/

#[cfg(test)]
mod tests {
    use super::*;
    use bio_seq::prelude::*;

    #[test]
    fn from_contig() {
        let seq: Seq<Dna> = dna!("AAAAAAAAAAAAAAAAATAAAGAAAAAAAAAAT");
        let debg: DeBruijn<Dna, 6, usize> = DeBruijn::from_seq(&seq);
        debg.dump();
        assert!(false);
    }
}
