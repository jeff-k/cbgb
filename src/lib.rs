#![feature(generic_const_exprs)]

use core::hash::Hash;
use std::cmp::Eq;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;

use bio_seq::codec::Codec;
use bio_seq::seq::iterators::KmerIter;
use bio_seq::{Kmer, Seq};

mod debruijn;

pub trait Monoid: Eq {
    fn zero() -> Self;
    fn one() -> Self;
    fn addm(&self, other: &Self) -> Self;
}

impl Monoid for usize {
    fn zero() -> Self {
        0
    }
    fn one() -> Self {
        1
    }
    fn addm(&self, other: &Self) -> Self {
        self + other
    }
}

trait GenomeGraph {
    type Edge;
    type Vertex;
    type Error;

    fn from_gfa() -> Result<Self, Self::Error>
    where
        Self: Sized;
    fn to_adj(self);
    fn subdawg(self);
    //    fn kmers(self) -> KmerIter;
    //    fn eulerian(self) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio_seq::codec::dna::Dna;
    use bio_seq::{dna, FromStr};

    #[test]
    fn fixme() {
        assert!(true);
    }
}
