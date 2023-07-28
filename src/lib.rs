use std::collections;
use std::ops::AddAssign;
use std::rc::Rc;

use hashbrown::{HashMap, HashSet};
use petgraph::visit::{
    GraphBase, GraphRef, IntoNeighbors, IntoNodeIdentifiers, NodeCount, NodeIndexable, Visitable,
};

use bio_seq::kmer::KmerIter;
use bio_seq::prelude::*;

pub struct GenomeGraph;

impl Default for GenomeGraph {
    fn default() -> Self {
        Self
    }
}

pub trait Edge: Default + Copy + PartialEq + Clone + From<u8> {}

impl Edge for u32 {}

#[derive(Clone)]
pub struct KmerIndex<E: Edge, const K: usize> {
    pub index: Vec<E>, // this could [E; { 4**K }] if const expression are allowed
    pub total: usize,
}

impl<E: Edge, const K: usize> GraphBase for KmerIndex<E, K> {
    type NodeId = Kmer<Dna, K>;
    type EdgeId = E;
}

/*
impl<E: Edge, const K: usize> GraphRef for KmerIndex<E, K> {
}
*/

#[derive(Clone)]
pub struct KmerTable<E: Edge, const K: usize> {
    pub index: HashMap<Kmer<Dna, K>, E>,
    pub total: usize,
}

impl<E: Edge, const K: usize> GraphBase for KmerTable<E, K> {
    type NodeId = Kmer<Dna, K>;
    type EdgeId = E;
}

/*
impl<'a, E: Edge, const K: usize> GraphRef for &'a KmerTable<E, K> {
}
*/

/*
impl<E: Edge, const K: usize> GraphRef for HashMap<Kmer<Dna, K>, E> {
}
*/

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

impl<'a, const K: usize> IntoNodeIdentifiers for &'a KmerIndex<u32, K> {
    type NodeIdentifiers = KmerIter<'a, Dna, K>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        unimplemented!()
    }
}

impl<'a, const K: usize> IntoNeighbors for &'a KmerIndex<u32, K> {
    type Neighbors = KmerIter<'a, Dna, K>;

    fn neighbors(self, kmer: Kmer<Dna, K>) -> Self::Neighbors {
        unimplemented!()
    }
}

impl<const K: usize> Visitable for KmerIndex<u32, K> {
    type Map = collections::HashSet<Kmer<Dna, K>>;
    fn visit_map(&self) -> Self::Map {
        unimplemented!()
    }

    fn reset_map(&self, map: &mut Self::Map) {
        unimplemented!()
    }
}

impl<E: Edge, const K: usize> Default for KmerIndex<E, K> {
    fn default() -> Self {
        KmerIndex {
            index: vec![E::default(); 1 << (K * Dna::WIDTH as usize)],
            total: 0,
        }
    }
}

impl<E: Edge + AddAssign<E>, const K: usize> Debruijn<K> for KmerIndex<E, K>
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

pub trait Debruijn<const K: usize>: GraphBase {
    fn add(&mut self, kmer: Kmer<Dna, K>);
    fn walk(&self, start: Kmer<Dna, K>);
    fn compress(&self) -> GenomeGraph;
    fn entropy(&self) -> f64;
    fn kld(&self, other: &Self) -> f64;
}

impl<E: Edge + AddAssign, const K: usize> Debruijn<K> for KmerTable<E, K>
where
    f64: From<E>,
{
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        let entry = self.index.entry(kmer).or_insert(E::default());
        entry.add_assign(E::from(1u8));
        self.total += 1;
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn compress(&self) -> GenomeGraph {
        unimplemented!()
    }

    fn entropy(&self) -> f64 {
        let mut h: f64 = 0.0;
        let mut t: f64 = 0.0;
        for count in self.index.values() {
            t += f64::from(*count);
        }
        for count in self.index.values() {
            let p: f64 = f64::from(*count) / t;
            if p > 0.0 {
                h += p * p.ln();
            }
        }
        h
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
