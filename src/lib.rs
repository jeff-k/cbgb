use std::collections::HashMap;
use std::ops::AddAssign;

use bio_seq::prelude::*;

pub struct GenomeGraph;

impl Default for GenomeGraph {
    fn default() -> Self {
        Self
    }
}

pub struct KmerSet;

impl Default for KmerSet {
    fn default() -> Self {
        Self
    }
}

pub struct KmerIndex<E, const K: usize> {
    pub index: Vec<E>,
    pub total: usize,
}

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
        let mut visited: KmerSet = KmerSet::default();
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
