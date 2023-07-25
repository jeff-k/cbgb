use std::ops::AddAssign;

use hashbrown::{HashMap, HashSet};
use petgraph::visit::{
    GraphBase, IntoNeighbors, IntoNodeIdentifiers, NodeCount, NodeIndexable, Visitable,
};

use bio_seq::prelude::*;

pub struct GenomeGraph;

impl Default for GenomeGraph {
    fn default() -> Self {
        Self
    }
}

pub struct KmerIndex<E, const K: usize> {
    pub index: Vec<E>,
    pub total: usize,
}

impl<E: Copy + PartialEq, const K: usize> GraphBase for KmerIndex<E, K> {
    type NodeId = Kmer<Dna, K>;
    type EdgeId = E;
}

impl<const K: usize> NodeCount for KmerIndex<u32, K> {
    fn node_count(&self) -> usize {
        self.index.len()
    }
}

impl<const K: usize> NodeIndexable for KmerIndex<u32, K> {
    fn node_bound(&self) -> usize {
        self.index.len()
    }

    fn to_index(&self, kmer: Kmer<Dna, K>) -> usize {
        usize::from(kmer)
    }

    fn from_index(&self, index: usize) -> Kmer<Dna, K> {
        unimplemented!()
    }
}

/*
impl<const K: usize> IntoNodeIdentifiers for KmerIndex<u32, K> {
    type NodeIdentifiers = Kmer<Dna, K>;

    fn node_identifiers(&self) -> Self::NodeIdentifiers {
        unimplemented!()
    }
}
*/

/*
impl<const K: usize> IntoNeighbors for KmerIndex<u32, K> {
    fn neighbors(&self, kmer: Kmer<Dna, K>) -> Self::Neighbors {
        unimplemented!()
    }
}
*/

/*
impl<const K: usize> Visitable for KmerIndex<u32, K> {
    type Map = HashMap<Kmer<Dna, K>, usize>;
    fn visit_map(&self) -> Self::Map {
        unimplemented!()
    }

    fn reset_map(&self, map: &mut Self::Map) {
        unimplemented!()
    }
}
*/

impl<E: Default + Copy, const K: usize> Default for KmerIndex<E, K> {
    fn default() -> Self {
        KmerIndex {
            index: vec![E::default(); 1 << (K * Dna::WIDTH as usize)],
            total: 0,
        }
    }
}

impl<E: Default + Copy + AddAssign<E> + From<u8>, const K: usize> Debruijn<K> for KmerIndex<E, K>
where
    f64: From<E>,
{
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        self.index[usize::from(kmer)].add_assign(E::from(1u8));
        self.total += 1;
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn compress(&self) -> GenomeGraph {
        //        let start: SeqSlice<Dna> = node.heads();
        let mut g: GenomeGraph = GenomeGraph::default();
        let mut visited: HashSet<Kmer<Dna, K>> = HashSet::new();
        g
    }

    fn entropy(&self) -> f64 {
        let mut h: f64 = 0.0;
        let t: f64 = self.total as f64;
        for i in 0..self.index.len() {
            let p: f64 = f64::from(self.index[i]) / t;
            if p > 0.0 {
                h += p * p.ln();
            }
        }
        h
    }

    fn kld(&self, other: &Self) -> f64 {
        let mut h: f64 = 0.0;
        let t: f64 = self.total as f64;
        for i in 0..self.index.len() {
            let p: f64 = f64::from(self.index[i]) / t;
            let q: f64 = f64::from(other.index[i]) / t;
            if p > 0.0 && q > 0.0 {
                h += q * p.ln();
            }
        }
        h
    }
}

pub trait Debruijn<const K: usize> {
    fn add(&mut self, kmer: Kmer<Dna, K>);
    fn walk(&self, start: Kmer<Dna, K>);
    fn compress(&self) -> GenomeGraph;
    fn entropy(&self) -> f64;
    fn kld(&self, other: &Self) -> f64;
}

impl<E: Default, const K: usize> Debruijn<K> for HashMap<Kmer<Dna, K>, E> {
    fn add(&mut self, _kmer: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn compress(&self) -> GenomeGraph {
        unimplemented!()
    }

    fn entropy(&self) -> f64 {
        unimplemented!()
    }

    fn kld(&self, _other: &Self) -> f64 {
        unimplemented!()
    }
}

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
