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

fn extend_condition<A: QuasiAlignment>(
    alignment: &A,
    r_pos: u32,
    q_pos: u32,
    forward: bool,
    k: u32,
) -> bool {
    if forward {
        alignment.r_end() + 1 == r_pos + k && alignment.q_end() + 1 == q_pos + k
    } else {
        alignment.r_start() == r_pos + 1 && alignment.q_end() + 1 == q_pos + k
    }
}

pub fn merge_segments<A: QuasiAlignment + Debug>(matches: Vec<Option<i32>>, k: u32) -> Vec<A> {
    let mut alignments: Vec<A> = Vec::new();

    // working values for the current segment
    let mut alignment: Option<A> = None;

    for (q_pos, pos) in matches.into_iter().enumerate() {
        let q_pos = q_pos as u32 + 1; // 1-indexed
        println!(
            "processing match q_pos: {:?}, pos: {:?}. Working alignment: {:?}",
            q_pos, pos, &alignment
        );

        match pos {
            Some(0) => {
                // Ambiguously matching kmer
                // You can decide what to do here based on your specific requirements
            }
            Some(r_pos) => {
                let forward = r_pos > 0;

                let r_pos = if forward {
                    r_pos.unsigned_abs()
                } else {
                    r_pos.unsigned_abs() - k
                };

                if let Some(a) = &mut alignment {
                    if extend_condition(a, r_pos, q_pos, forward, k) {
                        println!("\textend condition met");
                        // Extend the current alignment
                        if forward {
                            a.set_r_end(r_pos + k);
                        } else {
                            a.set_r_start(r_pos);
                        }
                        a.set_q_end(q_pos + k);
                    } else {
                        // Save the current alignment and start a new one
                        println!("\tpushing working alignment and starting new one");
                        alignments.push(a.clone());
                        alignment = Some(A::new(q_pos, q_pos + k, r_pos, r_pos + k, forward));
                    }
                } else {
                    println!("\tno initialised alignments, starting new one");
                    // Start a new alignment
                    alignment = Some(A::new(q_pos, q_pos + k, r_pos, r_pos + k, forward));
                }
            }
            None => {
                println!("\tNone for kmer match: resetting alignment");
                // Non-matching kmer
                if let Some(mut a) = alignment.take() {
                    // Save the current alignment if it exists
                    println!("\tpushing working alignment");
                    //   a.set_q_end(q_pos - 1);
                    alignments.push(a);
                }
                // Reset working alignment
                alignment = None;
            }
        }
    }

    // Save the last alignment if it exists
    if let Some(a) = alignment {
        alignments.push(a);
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
