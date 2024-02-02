use std::collections;
use std::ops::AddAssign;

use petgraph::visit::{
    GraphBase, IntoNeighbors, IntoNodeIdentifiers, NodeCount, NodeIndexable, Visitable,
};

use bio_seq::kmer::KmerIter;
use bio_seq::prelude::*;

use crate::{Debruijn, Edge, GenomeGraph, graph};

#[derive(Clone)]
pub struct KmerArray<E: Edge, const K: usize> {
    pub index: Vec<E>, // this could [E; { 4**K }] if const expression are allowed
    pub total: usize,
}

impl<E: Edge, const K: usize> GraphBase for KmerArray<E, K> {
    type NodeId = Kmer<Dna, K>;
    type EdgeId = E;
}

/*
impl<E: Edge, const K: usize> GraphRef for KmerArray<E, K> {
}
*/

impl<const K: usize> NodeCount for KmerArray<u32, K> {
    fn node_count(&self) -> usize {
        self.index.len()
    }
}

impl<const K: usize> NodeIndexable for KmerArray<u32, K> {
    fn node_bound(&self) -> usize {
        self.index.len()
    }

    fn to_index(&self, kmer: Kmer<Dna, K>) -> usize {
        usize::from(kmer)
    }

    fn from_index(&self, _index: usize) -> Kmer<Dna, K> {
        unimplemented!()
    }
}

impl<'a, const K: usize> IntoNodeIdentifiers for &'a KmerArray<u32, K> {
    type NodeIdentifiers = KmerIter<'a, Dna, K>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        unimplemented!()
    }
}

impl<'a, const K: usize> IntoNeighbors for &'a KmerArray<u32, K> {
    type Neighbors = KmerIter<'a, Dna, K>;

    fn neighbors(self, _kmer: Kmer<Dna, K>) -> Self::Neighbors {
        unimplemented!()
    }
}

impl<const K: usize> Visitable for KmerArray<u32, K> {
    type Map = collections::HashSet<Kmer<Dna, K>>; // would prefer hashbrown::HashSet
    fn visit_map(&self) -> Self::Map {
        unimplemented!()
    }

    fn reset_map(&self, _map: &mut Self::Map) {
        unimplemented!()
    }
}

impl<E: Edge, const K: usize> Default for KmerArray<E, K> {
    fn default() -> Self {
        KmerArray {
            index: vec![E::default(); 1 << (K * Dna::WIDTH as usize)],
            total: 0,
        }
    }
}

impl<E: Edge + AddAssign<E>, const K: usize> GenomeGraph for KmerArray<E, K> {
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        self.index[usize::from(kmer)].add_assign(E::from(1u8));
        self.total += 1;
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }
}

impl<E: Edge + AddAssign<E>, const K: usize> Debruijn<K> for KmerArray<E, K>
where
    f64: From<E>,
{
    fn compress(
        &self,
    ) -> graph::HashGraph
    {
        unimplemented!()
        //        let start: SeqSlice<Dna> = node.heads();
        //let g: GenomeGraph = GenomeGraph::default();
        //g
    }

    fn eulerian(&self) -> bool {
        unimplemented!()
    }
}

impl<E: Edge + AddAssign<E>, const K: usize> KmerArray<E, K> {
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
