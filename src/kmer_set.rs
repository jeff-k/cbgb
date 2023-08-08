use hashbrown::HashSet;

use bio_seq::prelude::*;

//use crate::{Debruijn, Edge, GenomeGraph};

#[derive(Clone)]
pub struct KmerSet<const K: usize> {
    pub index: HashSet<Kmer<Dna, K>>,
}

impl<const K: usize> Default for KmerSet<K> {
    fn default() -> Self {
        KmerSet {
            index: HashSet::new(),
        }
    }
}

impl<const K: usize> KmerSet<K> {
    pub fn contains(&self, kmer: Kmer<Dna, K>) -> bool {
        self.index.contains(&kmer)
    }

    pub fn add(&mut self, kmer: Kmer<Dna, K>) {
        self.index.insert(kmer);
    }

    pub fn len(&self) -> usize {
        self.index.len()
    }
}

/*
impl<const K: usize> Debruijn<K> for KmerSet<K>
{
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        self.index.insert(kmer);
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

    fn kld(&self, other: &Self) -> f64 {
        unimplemented!()
    }
}
*/
