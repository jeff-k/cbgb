mod kmer_array;
mod kmer_table;

pub use kmer_array::KmerArray;
pub use kmer_table::KmerTable;

use petgraph::visit::GraphBase;

//use bio_seq::kmer::KmerIter;
use bio_seq::prelude::*;

pub struct GenomeGraph;

impl Default for GenomeGraph {
    fn default() -> Self {
        Self
    }
}

pub trait Edge: Default + Copy + PartialEq + Clone + From<u8> {}

impl Edge for u32 {}
impl Edge for u16 {}
impl Edge for u8 {}

pub trait Debruijn<const K: usize>: GraphBase {
    fn add(&mut self, kmer: Kmer<Dna, K>);
    fn walk(&self, start: Kmer<Dna, K>);
    fn compress(&self) -> GenomeGraph;
    fn entropy(&self) -> f64;
    fn kld(&self, other: &Self) -> f64;
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
