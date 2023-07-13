use std::collections::HashMap;

use bio_seq::prelude::*;

pub trait Monoid {
    const M0: Self;
    const M1: Self;
    fn add(&mut self, rhs: Self);
}

impl Monoid for u16 {
    const M0: Self = 0u16;
    const M1: Self = 1u16;
    fn add(&mut self, rhs: Self) {
        *self += rhs
    }
}

impl Monoid for u32 {
    const M0: Self = 0u32;
    const M1: Self = 1u32;
    fn add(&mut self, rhs: Self) {
        *self += rhs
    }
}

pub struct GenomeGraph;

pub struct KmerIndex<E, const K: usize> {
    pub index: Vec<E>,
    pub total: E,
}

impl<E: Monoid + Copy, const K: usize> Default for KmerIndex<E, K> {
    fn default() -> Self {
        KmerIndex {
            index: vec![E::M0; 1 << (K * Dna::WIDTH as usize)],
            total: E::M0,
        }
    }
}

impl<E: Monoid + Copy, const K: usize> Debruijn<K> for KmerIndex<E, K>
where
    f32: From<E>,
{
    fn add(&mut self, kmer: Kmer<Dna, K>) {
        self.index[usize::from(kmer)].add(E::M1);
        self.total.add(E::M1);
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn compress(&self) -> GenomeGraph {
        unimplemented!()
    }

    fn entropy(&self) -> f32 {
        let mut h: f32 = 0.0;
        for i in 0..self.index.len() {
            let p: f32 = f32::from(self.index[i]) / f32::from(self.total);
            if p > 0.0 {
                h += p * p.ln();
            }
        }
        h
    }

    fn kld(&self, other: &Self) -> f32 {
        let mut h: f32 = 0.0;
        for i in 0..self.index.len() {
            let p: f32 = f32::from(self.index[i]) / f32::from(self.total);
            let q: f32 = f32::from(other.index[i]) / f32::from(self.total);
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
    fn entropy(&self) -> f32;
    fn kld(&self, other: &Self) -> f32;
}

impl<E: Monoid + Copy, const K: usize> Debruijn<K> for HashMap<Kmer<Dna, K>, E> {
    fn add(&mut self, _kmer: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn walk(&self, _start: Kmer<Dna, K>) {
        unimplemented!()
    }

    fn compress(&self) -> GenomeGraph {
        unimplemented!()
    }

    fn entropy(&self) -> f32 {
        unimplemented!()
    }

    fn kld(&self, _other: &Self) -> f32 {
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
