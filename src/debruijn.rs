#![feature(generic_const_exprs)]

use core::hash::Hash;
use std::cmp::Eq;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;

use bio_seq::codec::Codec;
use bio_seq::seq::iterators::KmerIter;
use bio_seq::{Kmer, Seq};

use petgraph::visit::{GraphBase, GraphRef, IntoNeighbors};

use crate::{GenomeGraph, Monoid};

pub struct DeBruijn<A: Codec + Eq, const N: usize, Edge: Hash + Monoid> {
    map: HashMap<Kmer<A, { N + 1 }>, Edge>,
}

/*
pub struct DenseDeBruijn<A: Codec, const N: usize, Edge: Hash + Monoid> {
    map: Vec<Edge>,
}
*/

impl<A: Codec + Eq + Display, const N: usize, Edge: Monoid + Hash + Display> DeBruijn<A, N, Edge>
where
    [(); N + 1]:,
{
    /*
    fn compress(self) -> GenomeGraph<A, Edge> {
        GenomeGraph {
            edges: HashMap::new(),
            vertices: HashSet::new(),
        }
    }
    */
    fn from_kmers(iter: KmerIter<A, { N + 1 }>) -> Self {
        let mut kmers = HashMap::new();
        for kmer in iter {
            let mut r = kmers.entry(kmer).or_insert(Edge::zero());
            *r = r.addm(&Edge::one());
        }
        DeBruijn { map: kmers }
    }
    fn push_kmer(mut self, kmer: Kmer<A, { N + 1 }>) {
        self.map.insert(kmer, Edge::zero());
    }
    fn dump(self) {
        let mut i = 0;
        for (kmer, edge) in self.map {
            println!("{i}\t{kmer}\t{edge}");
            i += 1;
        }
    }
}

impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphBase for DeBruijn<A, N, Edge>
where
    [(); N + 1]:,
{
    type EdgeId = Edge;
    type NodeId = Kmer<A, N>;
}

/*
impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphRef for &DeBruijn<A, N, Edge>
where [(); N+1]:,{}
*/

impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> IntoNeighbors
    for &DeBruijn<A, N, Edge>
where
    [(); N + 1]:,
{
    type Neighbors = KmerIter<A, N>;

    fn neighbors(self) -> Self::Neighbors {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio_seq::codec::dna::Dna;
    use bio_seq::{dna, FromStr};

    #[test]
    fn from_contig() {
        let seq: Seq<Dna> = dna!("AAAAAAAAAAAAAAAAATAAAGAAAAAAAAAAT");
        let debg: DeBruijn<Dna, 6, usize> = DeBruijn::from_kmers(seq.kmers());
        debg.dump();
        assert!(false);
    }
}
