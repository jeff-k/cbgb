use core::hash::Hash;
use std::cmp::Eq;
use std::collections::HashMap;
//use std::collections::HashSet;
use std::fmt::{Debug, Display};

use bio_seq::codec::Codec;
//use bio_seq::seq::iterators::KmerIter;
use bio_seq::{Kmer, Seq};

use petgraph::visit::{GraphBase, IntoNeighbors};

use crate::Monoid;

pub struct DeBruijn<A: Codec + Eq, const N: usize, Edge: Hash + Monoid> {
    map: HashMap<Kmer<A, N>, (A, Edge)>,
}

/*
pub struct DenseDeBruijn<A: Codec, const N: usize, Edge: Hash + Monoid> {
    map: Vec<Edge>,
}
*/

impl<A: Codec + Eq + Display + Debug, const N: usize, Edge: Monoid + Hash + Display>
    DeBruijn<A, N, Edge>
{
    /*
    fn compress(self) -> GenomeGraph<A, Edge> {
        GenomeGraph {
            edges: HashMap::new(),
            vertices: HashSet::new(),
        }
    }
    */
    fn from_seq(seq: Seq<A>) -> Self {
        if seq.len() < N {
            panic!("Error: seq too short");
        }
        let mut kmers = HashMap::new();
        for subseq in seq.windows(N + 1) {
            let kin: Kmer<A, N> = subseq[..N].into();
            let kout: Kmer<A, N> = subseq[1..].into();
            //let mut r = kmers.entry(kmer).or_insert((out, Edge::zero()));
            //*r = r.addm(&Edge::one());
        }
        DeBruijn { map: kmers }
    }
    /*
    fn push_kmer(mut self, kmer: Kmer<A, N>) {
        self.map.insert(kmer, (Edge::zero());
    }
    */
    fn dump(self) {
        let mut i = 0;
        for (kmer, edge) in self.map {
            println!("{i}\t{kmer}\t{edge}");
            i += 1;
        }
    }
}

impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphBase for DeBruijn<A, N, Edge> {
    type EdgeId = Edge;
    type NodeId = Kmer<A, N>;
}

/*
impl<A: Codec + Eq, const N: usize, Edge: Monoid + Hash + Copy> GraphRef for &DeBruijn<A, N, Edge>
where [(); N+1]:,{}
*/

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
