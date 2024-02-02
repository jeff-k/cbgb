use hashbrown::{hash_map, HashMap};
use petgraph::visit::{GraphBase, IntoEdgeReferences};

use bio_seq::prelude::*;

use crate::{Edge, GenomeGraph, Path};

#[derive(Clone)]
pub struct HashGraph {
    pub index: Vec<Seq<Dna>>,
    pub graph: HashMap<usize, usize>,
}

impl<'a> IntoEdgeReferences for &'a HashGraph {
    type EdgeRef = usize;
    type EdgeReferences = hash_map::Values<'a, usize, usize>;

    fn edge_references(self) -> Self::EdgeReferences {
        self.edges.values()
    }
}

impl GraphBase for HashGraph {
    type NodeId = usize;
    type EdgeId = usize;
}

impl Default for HashGraph {
    fn default() -> Self {
        unimplemented!()
    }
}

impl GenomeGraph for HashGraph {
    fn add(&mut self, seq: &SeqSlice<Dna>) {
        // push copy of seq into index Vec, and add it's index to
        // the graph hashmap
        unimplemented!()
    }

    fn walk<'a>(&'a self, start: Self::NodeId) -> Path<'a, Self> {
        unimplemented!()
    }
}
