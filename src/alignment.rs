use bio_seq::prelude::{Dna, SeqSlice};
use core::fmt::Debug;

pub trait QuasiAlignment: Clone {
    fn new(q_start: u32, q_end: u32, r_start: u32, r_end: u32, forward: bool) -> Self;
    fn q_start(&self) -> u32;
    fn q_end(&self) -> u32;
    fn r_start(&self) -> u32;
    fn r_end(&self) -> u32;
    fn forward(&self) -> bool;
    fn set_q_start(&mut self, q_start: u32);
    fn set_q_end(&mut self, q_end: u32);
    fn set_r_start(&mut self, r_start: u32);
    fn set_r_end(&mut self, r_end: u32);
    fn set_forward(&mut self, forward: bool);
}

#[derive(Clone, Debug)]
pub struct Alignment {
    q_start: u32,
    q_end: u32,
    r_start: u32,
    r_end: u32,
    forward: bool,
}

impl QuasiAlignment for Alignment {
    fn new(q_start: u32, q_end: u32, r_start: u32, r_end: u32, forward: bool) -> Self {
        Alignment {
            q_start,
            q_end,
            r_start,
            r_end,
            forward,
        }
    }

    fn q_start(&self) -> u32 {
        self.q_start
    }

    fn q_end(&self) -> u32 {
        self.q_end
    }

    fn r_start(&self) -> u32 {
        self.r_start
    }

    fn r_end(&self) -> u32 {
        self.r_end
    }

    fn forward(&self) -> bool {
        self.forward
    }

    fn set_q_start(&mut self, q_start: u32) {
        self.q_start = q_start;
    }

    fn set_q_end(&mut self, q_end: u32) {
        self.q_end = q_end;
    }

    fn set_r_start(&mut self, r_start: u32) {
        self.r_start = r_start;
    }

    fn set_r_end(&mut self, r_end: u32) {
        self.r_end = r_end;
    }

    fn set_forward(&mut self, forward: bool) {
        self.forward = forward;
    }
}

pub trait QuasiAlign<A: QuasiAlignment> {
    fn quasi_align(&self, seq: &SeqSlice<Dna>) -> Vec<A>;
}

pub fn mergable<A: QuasiAlignment>(working: &A, current: &A) -> bool {
    working.forward() == current.forward()
        && working.r_start() + 1 == current.r_start()
        && working.q_end() + 1 == current.q_start()
}

pub fn merge<A: QuasiAlignment>(working: &mut A, current: &A) {
    working.set_q_end(current.q_end());
    working.set_r_end(current.r_end());
}

pub fn merge_segments<A: QuasiAlignment + Debug>(matches: Vec<Option<i32>>, k: u32) -> Vec<A> {
    let mut alignments: Vec<A> = Vec::new();
    let mut alignment: Option<A> = None; // the working alignment

    for (q_start, &mapping) in matches.iter().enumerate() {
        let q_start: u32 = q_start as u32 + 1; // 1-indexed
        let q_end: u32 = q_start + k;

        let segment: Option<A> = match mapping {
            Some(0) => {
                // if this kmer matches the reference ambiguously
                None
            }
            Some(r_pos) => {
                // if this kmer matches the reference uniquely
                let forward: bool = r_pos > 0;

                let r_start: u32 = if forward {
                    r_pos as u32
                } else {
                    r_pos.unsigned_abs() - k
                };
                let r_end: u32 = r_start + k;
                Some(A::new(q_start, q_end, r_start, r_end, forward))
            }
            None => {
                // if this kmer does not match the reference
                None
            }
        };

        println!("{:?}", segment);

        match (&mut alignment, &segment) {
            (None, None) => {
                println!("\t(no working alignment or segment)");
            }
            (Some(working_alignment), None) => {
                // encountering a run of None, push alignment
                println!("\tworking alignment but segment None");
                alignments.push(working_alignment.clone());
                alignment = None;
            }
            (Some(working_alignment), Some(current_segment)) => {
                // decide whether to merge current segment into working alignment
                if mergable(working_alignment, current_segment) {
                    println!("\tmergable segment");
                    merge(working_alignment, current_segment);
                } else {
                    println!("\tnon-mergeable segment. new alignment");
                    alignments.push(working_alignment.clone());
                    alignment = Some(current_segment.clone());
                }
            }
            (None, Some(current_segment)) => {
                // initialize new working alignment from this segment
                println!("\tNo working alignment, new segment");
                alignment = Some(current_segment.clone());
            }
        }
    }

    // if there is a remaining working alignment, push it
    if let Some(remaining) = alignment {
        alignments.push(remaining);
    }
    alignments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let kmer_map = vec![Some(1), Some(2), Some(3), Some(4)];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        assert_eq!(alignments.len(), 1);
        let a = &alignments[0];
        assert_eq!(a.q_start, 1);
        assert_eq!(a.q_end, 4);
        assert_eq!(a.r_start, 1);
        assert_eq!(a.r_end, 4);
        assert!(a.forward);
    }

    #[test]
    fn test_non_matching_kmers() {
        let kmer_map = vec![None, None, None];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        assert_eq!(alignments.len(), 0);
    }

    #[test]
    fn test_ambiguous_kmers() {
        let kmer_map = vec![Some(0), Some(0), Some(0)];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        // Assuming you don't create alignments for ambiguous kmers
        assert_eq!(alignments.len(), 0);
    }

    #[test]
    fn test_multiple_alignments() {
        let kmer_map = vec![Some(1), Some(2), None, Some(3), Some(4)];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        assert_eq!(alignments.len(), 2);
        let a = &alignments[0];
        assert_eq!(a.q_start, 1);
        assert_eq!(a.q_end, 2);
        assert_eq!(a.r_start, 1);
        assert_eq!(a.r_end, 2);
        assert!(a.forward);
        let b = &alignments[1];
        assert_eq!(b.q_start, 4);
        assert_eq!(b.q_end, 5);
        assert_eq!(b.r_start, 3);
        assert_eq!(b.r_end, 4);
        assert!(b.forward);
    }

    #[test]
    fn test_reverse_alignment() {
        let kmer_map = vec![Some(-4), Some(-3), Some(-2), Some(-1)];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        assert_eq!(alignments.len(), 1);
        let a = &alignments[0];
        assert_eq!(a.q_start, 1);
        assert_eq!(a.q_end, 4);
        assert_eq!(a.r_start, 1);
        assert_eq!(a.r_end, 4);
        assert!(a.forward == false);
    }

    #[test]
    fn test_multiple_reverse_alignments() {
        let kmer_map = vec![Some(-4), Some(-3), None, Some(-2), Some(-1)];
        let alignments: Vec<Alignment> = merge_segments(kmer_map, 0);
        assert_eq!(alignments.len(), 2);
        let a = &alignments[0];
        assert_eq!(a.q_start, 1);
        assert_eq!(a.q_end, 2);
        assert_eq!(a.r_start, 3);
        assert_eq!(a.r_end, 4);
        assert!(a.forward == false);
        let b = &alignments[1];
        assert_eq!(b.q_start, 4);
        assert_eq!(b.q_end, 5);
        assert_eq!(b.r_start, 1);
        assert_eq!(b.r_end, 2);
        assert!(b.forward == false);
    }

    #[test]
    fn test_self_alignment_with_k_5() {
        let kmer_map = vec![
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            Some(7),
            Some(8),
            Some(9),
            Some(10),
            Some(11),
            Some(12),
            Some(13),
            Some(14),
            Some(15),
            Some(16),
            Some(17),
            Some(18),
            Some(19),
            Some(20),
        ];
        let k = 4; // K-mer length is 5, so k = 5 - 1 = 4
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        assert_eq!(alignments.len(), 1);
        let a = &alignments[0];
        assert_eq!(a.q_start(), 1);
        assert_eq!(a.q_end(), 24);
        assert_eq!(a.r_start(), 1);
        assert_eq!(a.r_end(), 24);
        assert!(a.forward());
    }

    #[test]
    fn test_self_reverse_alignment_with_k_5() {
        let kmer_map = vec![
            Some(-24),
            Some(-23),
            Some(-22),
            Some(-21),
            Some(-20),
            Some(-19),
            Some(-18),
            Some(-17),
            Some(-16),
            Some(-15),
            Some(-14),
            Some(-13),
            Some(-12),
            Some(-11),
            Some(-10),
            Some(-9),
            Some(-8),
            Some(-7),
            Some(-6),
            Some(-5),
        ];
        let k = 4; // K-mer length is 5, so k = 5 - 1 = 4
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 1);
        let a = &alignments[0];
        assert_eq!(a.q_start(), 1);
        assert_eq!(a.q_end(), 24);
        assert_eq!(a.r_start(), 1);
        assert_eq!(a.r_end(), 24);
        assert!(!a.forward());
    }

    #[test]
    fn test_repeats_forward() {
        let kmer_map = vec![
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
        ];
        let k = 4; // K-mer length is 5, so k = 5 - 1 = 4
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        assert_eq!(alignments.len(), 3);

        let a = &alignments[0];
        assert_eq!(a.q_start(), 1);
        assert_eq!(a.q_end(), 6);
        assert_eq!(a.r_start(), 1);
        assert_eq!(a.r_end(), 6);
        assert!(a.forward());
        let b = &alignments[1];
        assert_eq!(b.q_start(), 7);
        assert_eq!(b.q_end(), 12);
        assert_eq!(b.r_start(), 1);
        assert_eq!(b.r_end(), 6);
        assert!(b.forward());
        let c = &alignments[2];
        assert_eq!(c.q_start(), 13);
        assert_eq!(c.q_end(), 18);
        assert_eq!(c.r_start(), 1);
        assert_eq!(c.r_end(), 6);
        assert!(c.forward());
    }

    #[test]
    fn test_unmapped_ends() {
        let kmer_map = vec![
            None,
            None,
            None,
            None,
            None,
            None,
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            Some(7),
            Some(8),
            None,
            None,
            None,
            None,
            None,
            None,
        ];
        let k = 4; // K-mer length is 5, so k = 5 - 1 = 4
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        assert_eq!(alignments.len(), 1);

        let a = &alignments[0];
        assert_eq!(a.q_start(), 7);
        assert_eq!(a.q_end(), 14);
        assert_eq!(a.r_start(), 1);
        assert_eq!(a.r_end(), 8);
        assert!(a.forward());
    }

    #[test]
    fn test_disjoint_segments_with_k_5() {
        let kmer_map = vec![
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5), // Forward segment
            None,
            None,
            None,
            None,
            None, // Gap
            Some(11),
            Some(12),
            Some(13),
            Some(14),
            Some(15),
            Some(16),
            Some(17),
            Some(18),
            Some(19),
            Some(20),
            Some(21), // Forward segment
            Some(-24),
            Some(-23),
            Some(-22),
            Some(-21),
            Some(-20), // Reverse segment
        ];
        let k = 4; // K-mer length is 5, so k = 5 - 1 = 4
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        assert_eq!(alignments.len(), 3); // Two disjoint segments in forward and one in reverse

        let a = &alignments[0];
        assert_eq!(a.q_start(), 1);
        assert_eq!(a.q_end(), 5);
        assert_eq!(a.r_start(), 1);
        assert_eq!(a.r_end(), 5);
        assert!(a.forward());

        let b = &alignments[1];
        assert_eq!(b.q_start(), 11);
        assert_eq!(b.q_end(), 21);
        assert_eq!(b.r_start(), 11);
        assert_eq!(b.r_end(), 21);
        assert!(b.forward());

        let c = &alignments[2];
        assert_eq!(c.q_start(), 22);
        assert_eq!(c.q_end(), 26);
        assert_eq!(c.r_start(), 20);
        assert_eq!(c.r_end(), 24);
        assert!(!c.forward());
    }
}
