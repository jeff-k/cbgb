use hashbrown::HashMap;

use bio_seq::prelude::*;

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
            index: index,
            len: len as u32,
        }
    }

    pub fn map(&self, seq: &SeqSlice<Dna>) -> Vec<Option<i32>> {
        if seq.len() < K {
            panic!()
        }

        let mut mapping: Vec<Option<i32>> = Vec::new();
        for kmer in seq.kmers() {
            mapping.push(self.index.get(&kmer).copied());
        }
        mapping
    }
}
