use std::cmp::Eq;

mod debruijn;

pub trait Monoid: Eq {
    const ZERO: Self;
    const ONE: Self;
    fn addm(&self, other: &Self) -> Self;
}

impl Monoid for () {
    const ZERO: () = ();
    const ONE: () = ();

    fn addm(&self, _other: &Self) -> Self {}
}

impl Monoid for usize {
    const ZERO: usize = 0;
    const ONE: usize = 1;

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
