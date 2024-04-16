pub mod alignment;
//mod graph;
//mod kmer_array;
mod kmer_map;
mod kmer_set;
//mod kmer_table;

pub use alignment::{Alignment, QuasiAlign};
//pub use kmer_array::KmerArray;
pub use kmer_map::KmerMap;
pub use kmer_set::KmerSet;
//pub use kmer_set::KmerSet;
//pub use kmer_table::KmerTable;

//use petgraph::visit::{GraphBase, IntoEdgeReferences};

use bio_seq::prelude::*;

pub trait Edge: Default + Copy + PartialEq + Clone + From<bool> {}

impl Edge for u32 {}
impl Edge for u16 {}
impl Edge for u8 {}
impl Edge for usize {}
impl Edge for bool {}

//pub struct Path<'a, G: GenomeGraph> {
//    graph: &'a G,
//    edges: Vec<G::EdgeRef>,
//}

//pub trait GenomeGraph: GraphBase + IntoEdgeReferences {
//    fn add(&mut self, seq: &SeqSlice<Dna>);
//    fn walk(&self, start: Self::NodeId) -> Path<Self>;
//}

//pub trait Debruijn<const K: usize>: GenomeGraph {
//    fn compress(&self) -> graph::HashGraph;
//    fn eulerian(&self) -> bool;
//}

pub trait RankSelect {
    fn rank(reference: &SeqSlice<Dna>, query: &SeqSlice<Dna>) -> usize;
    fn select(reference: &SeqSlice<Dna>, query: &SeqSlice<Dna>, rank: usize) -> usize;
}
