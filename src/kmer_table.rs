use std::ops::AddAssign;

use hashbrown::HashMap;
use petgraph::visit::GraphBase;

//use bio_seq::kmer::KmerIter;
use bio_seq::prelude::*;

use crate::{graph, Debruijn, Edge, GenomeGraph};

#[derive(Clone)]
pub struct KmerTable<E: Edge, const K: usize> {
    pub index: HashMap<Kmer<Dna, K>, E>,
    pub total: usize,
}

impl<E: Edge, const K: usize> GraphBase for KmerTable<E, K> {
    type NodeId = Kmer<Dna, K>;
    type EdgeId = E;
}

impl<E: Edge, const K: usize> Default for KmerTable<E, K> {
    fn default() -> Self {
        KmerTable {
            index: HashMap::new(),
            total: 0,
        }
    }
}

/*
impl<'a, E: Edge, const K: usize> GraphRef for &'a KmerTable<E, K> {
}
*/

impl<E: Edge + AddAssign, const K: usize> GenomeGraph for KmerTable<E, K>
where
    f64: From<E>,
{
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        let entry = self.index.entry(kmer).or_default();
        entry.add_assign(E::from(1u8));
        self.total += 1;
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }
}

impl<E: Edge + AddAssign, const K: usize> Debruijn<K> for KmerTable<E, K>
where
    f64: From<E>,
{
    fn compress(&self) -> graph::HashGraph {
        unimplemented!()
    }

    fn eulerian(&self) -> bool {
        unimplemented!()
    }
}

impl<E: Edge + AddAssign, const K: usize> KmerTable<E, K>
where
    f64: From<E>,
{
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

    fn kld(&self, other: &Self) -> f64 {
        let mut h: f64 = 0.0;
        let t_p: f64 = self.total as f64;
        let t_q: f64 = other.total as f64;
        if t_p == 0.0 || t_q == 0.0 {
            return 0.0;
        }

        for (kmer, c_p) in self.index.iter() {
            if let Some(c_q) = other.index.get(kmer) {
                let p: f64 = f64::from(*c_p) / t_p;
                let q: f64 = f64::from(*c_q) / t_q;
                h += q * p.ln();
            }
        }
        h
    }
}
