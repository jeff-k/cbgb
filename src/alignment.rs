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

#[derive(Clone, Debug, PartialEq)]
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
    if working.forward() == current.forward() {
        if working.forward() {
            working.r_end() + 1 == current.r_end() // || working.r_end() == current.r_end())
                && (working.q_end() > current.q_start())
        } else {
            working.r_start() == current.r_start() + 1 // || working.r_start() == current.r_start())
                && (working.q_end() > current.q_start())
        }
    } else {
        false
    }
}

pub fn merge<A: QuasiAlignment>(working: &mut A, current: &A) {
    working.set_q_end(current.q_end());
    if working.forward() {
        working.set_r_end(current.r_end());
    } else {
        working.set_r_start(current.r_start());
    }
}

pub fn merge_segments<A: QuasiAlignment + Debug>(matches: Vec<Option<i32>>, k: u32) -> Vec<A> {
    let mut alignments: Vec<A> = Vec::new();
    let mut alignment: Option<A> = None; // the working alignment
    let mut last_segment_end: u32 = 0;

    let mut match_enumerator = matches.iter().enumerate();
    for (q_start, &mapping) in match_enumerator {
        let q_start: u32 = q_start as u32;
        let q_end: u32 = q_start + k;

        if q_start < last_segment_end {
            println!(
                "\t\tskipping {}-{}, {:?} (last q_end: {})",
                q_start, q_end, mapping, last_segment_end
            );
            continue;
        }

        let segment: Option<A> = match mapping {
            Some(0) => {
                // if this kmer matches the reference ambiguously
                unimplemented!()
            }
            Some(r_pos) => {
                // if this kmer matches the reference uniquely
                let forward: bool = r_pos > 0;

                let (r_start, r_end): (u32, u32) = if forward {
                    let r_pos: u32 = r_pos as u32 - 1;
                    (r_pos, r_pos + k)
                } else {
                    let r_pos: u32 = r_pos.unsigned_abs();
                    (r_pos - k, r_pos)
                };
                Some(A::new(q_start, q_end, r_start, r_end, forward))
            }
            None => {
                // if this kmer does not match the reference
                None
            }
        };

        println!("working alignment: {:?}", alignment);
        println!("new segment:       {:?}", segment);

        match (&mut alignment, &segment) {
            (None, None) => {
                println!("\t(no working alignment or segment)");
            }
            (Some(working_alignment), None) => {
                // encountering a run of None, push alignment
                println!("\tworking alignment but segment None");
                last_segment_end = working_alignment.q_end();
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
                    last_segment_end = working_alignment.q_end();
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
        println!("\tquery exhausted, pushing remaining working segment");
        alignments.push(remaining);
    }
    alignments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_alignment() {
        // Testing an alignment scenario where the query is a substring of the reference
        let kmer_map = vec![Some(7), Some(8), Some(9), Some(10), Some(11)];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 1);
        let a = &alignments[0];
        assert_eq!(
            *a,
            Alignment {
                q_start: 0,
                q_end: 9,
                r_start: 6,
                r_end: 15,
                forward: true,
            }
        );
    }

    #[test]
    fn test_multiple_alignments() {
        // Testing a scenario where the query has multiple regions that align to the reference
        // Original Reference: "ACGTGACGGTCGTACCACCAAAGT", Query: "ACGTGACGGTACGTGACGGT"
        // aligner.map_seq(query): [1, 2, 3, 4, 5, 6, None, -16, -15, None, 1, 2, 3, 4, 5, 6]
        let kmer_map = vec![
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            None,
            // this example generates a breakpoint artefact
            None, //Some(-16),
            None, //Some(-15),
            None,
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
        ];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        //assert_eq!(alignments.len(), 2);
        let a = &alignments[0];
        assert_eq!(
            *a,
            Alignment {
                q_start: 0,
                q_end: 10,
                r_start: 0,
                r_end: 10,
                forward: true
            }
        );
        let b = &alignments[1];
        assert_eq!(
            *b,
            Alignment {
                q_start: 10,
                q_end: 20,
                r_start: 0,
                r_end: 10,
                forward: true
            }
        );
    }

    #[test]
    fn test_repeats_in_sequence() {
        // Testing a scenario where the query has repeats that align to different parts of the reference
        // Original Reference: "ACGTGACGGTCGTACCACCAAAGT", Query: "ACCGTCACGTACGTGACGGT"
        // aligner.map_seq(query): [-10, -9, -8, -7, -6, -5, None, 11, -15, None, 1, 2, 3, 4, 5, 6]
        let kmer_map = vec![
            Some(-10),
            Some(-9),
            Some(-8),
            Some(-7),
            Some(-6),
            Some(-5),
            None,
            // this example generates a breakpoint artefact
            Some(11),
            Some(-15),
            None,
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
        ];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 2);
        let a = &alignments[0];
        assert_eq!(
            *a,
            Alignment {
                q_start: 0,
                q_end: 10,
                r_start: 0,
                r_end: 10,
                forward: false
            }
        );
        let b = &alignments[1];
        assert_eq!(
            *b,
            Alignment {
                q_start: 10,
                q_end: 20,
                r_start: 0,
                r_end: 10,
                forward: true
            }
        );
    }

    #[test]
    fn test_basic_alignment() {
        // Testing a basic alignment scenario

        let kmer_map = vec![Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 1);

        let a0 = &alignments[0];
        assert_eq!(
            *a0,
            Alignment {
                q_start: 0,
                q_end: 10,
                r_start: 0,
                r_end: 10,
                forward: true
            }
        );
    }

    #[test]
    fn test_basic_reverse_alignment() {
        // Testing a basic alignment scenario

        let kmer_map = vec![Some(-10), Some(-9), Some(-8), Some(-7), Some(-6), Some(-5)];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 1);

        let a0 = &alignments[0];
        assert_eq!(
            *a0,
            Alignment {
                q_start: 0,
                q_end: 10,
                r_start: 0,
                r_end: 10,
                forward: false
            }
        );
    }

    #[test]
    fn test_non_matching_kmers() {
        // Testing a scenario with non-matching k-mers
        // Original Reference: ACGTGACGGTCGTACCACCAAAGT
        // Query: ACGTNACGG
        // k: 5

        let kmer_map = vec![None, None, None, None, None];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 0);
    }

    #[test]
    fn test_ambiguous_kmers() {
        // Testing a scenario with ambiguous k-mers
        // note that the middle `G` should have matches belong to two
        // different segments.
        // Original Reference: ACGTGACGGTCGTACCACCAAAGT
        // Query:              ACGTG       TACCAC
        // k: 5

        let kmer_map = vec![Some(1), None, None, None, Some(12), Some(13), Some(14)];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        assert_eq!(alignments.len(), 2);

        let a0 = &alignments[0];
        assert_eq!(
            *a0,
            Alignment {
                q_start: 0,
                q_end: 5,
                r_start: 0,
                r_end: 5,
                forward: true
            }
        );

        let a1 = &alignments[1];
        assert_eq!(
            *a1,
            Alignment {
                q_start: 5,
                q_end: 11,
                r_start: 12,
                r_end: 18,
                forward: true
            }
        );
    }

    #[test]
    fn test_interrupted_alternating_direction() {
        // ..
        // Original Reference: ACGTGACGGTCGTACCACCAAAGT
        // Query:                      GGTACGAC
        //                               NGTACCACC
        // k: 5

        let kmer_map = vec![
            Some(-16),
            Some(-15),
            Some(-14),
            Some(-13),
            None,
            None,
            None,
            None,
            None,
            Some(12),
            Some(13),
            Some(14),
            Some(15),
        ];
        let k = 5; // K-mer length is 5
        let alignments: Vec<Alignment> = merge_segments(kmer_map, k);
        println!("{:?}", alignments);
        assert_eq!(alignments.len(), 2);
        let a0 = &alignments[0];
        assert_eq!(
            *a0,
            Alignment {
                q_start: 0,
                q_end: 8,
                r_start: 8,
                r_end: 16,
                forward: false
            }
        );

        let a1 = &alignments[1];
        assert_eq!(
            *a1,
            Alignment {
                q_start: 9,
                q_end: 17,
                r_start: 11,
                r_end: 19,
                forward: true
            }
        );
    }
}
