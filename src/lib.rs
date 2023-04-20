use core::hash::Hash;
use std::cmp::Eq;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;

mod debruijn;

pub trait Monoid: Eq {
    const zero: Self;
    const one: Self;
    fn addm(&self, other: &Self) -> Self;
}

impl Monoid for () {
    const zero: () = ();
    const one: () = ();

    fn addm(&self, other: &Self) -> Self {
        ()
    }
}

impl Monoid for usize {
    const zero: usize = 0;
    const one: usize = 1;

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

    #[test]
    fn fixme() {
        assert!(true);
    }
}
