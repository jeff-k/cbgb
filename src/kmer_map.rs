use hashbrown::HashMap;

use bio_seq::prelude::*;

use crate::alignment;

#[derive(Clone)]
pub struct KmerMap<const K: usize> {
    pub index: HashMap<Kmer<Dna, K>, i32>,
    pub len: u32,
}

impl<const K: usize> KmerMap<K> {
    pub fn new(seq: &SeqSlice<Dna>) -> Self {
        let mut index: HashMap<Kmer<Dna, K>, i32> = HashMap::new();
        let len: i32 = seq.len() as i32;

        for (pos, kmer) in seq.kmers().enumerate() {
            if let Some(v) = index.get_mut(&kmer) {
                *v = 0;
            } else {
                index.insert(kmer, pos as i32 + 1);
            }
        }
        for (pos, kmer) in seq.revcomp().kmers().enumerate() {
            if let Some(v) = index.get_mut(&kmer) {
                *v = 0;
            } else {
                index.insert(kmer, (pos as i32 - 1) - len);
            }
        }

        KmerMap {
            index,
            len: len as u32,
        }
    }

    pub fn match_kmers(&self, seq: &SeqSlice<Dna>) -> Vec<Option<i32>> {
        let mut mapping: Vec<Option<i32>> = Vec::new();

        if seq.len() < K {
            // this may better be an exception
            return mapping;
        }

        for kmer in seq.kmers() {
            mapping.push(self.index.get(&kmer).copied());
        }
        mapping
    }
}

impl<const K: usize, A: alignment::QuasiAlignment> alignment::QuasiAlign<A> for KmerMap<K> {
    fn quasi_align(&self, seq: &SeqSlice<Dna>) -> Vec<A> {
        let v: Vec<Option<i32>> = self.match_kmers(seq);
        alignment::merge_segments(v, K)
    }
}
