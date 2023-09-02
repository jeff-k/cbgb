use bio_seq::prelude::{Dna, SeqSlice};

#[derive(Clone, Debug)]
pub struct Alignment {
    q_start: u32,
    q_end: u32,
    r_start: u32,
    r_end: u32,
    forward: bool,
}

pub trait QuasiAlign {
    fn quasi_align(&self, seq: &SeqSlice<Dna>) -> Vec<Alignment>;
}

pub fn merge_segments(matches: Vec<Option<i32>>, _k: usize) -> Vec<Alignment> {
    let mut alignments: Vec<Alignment> = Vec::new();

    // working values for the current segment
    let mut alignment: Option<Alignment> = None;

    for (q_pos, pos) in matches.into_iter().enumerate() {
        let q_pos = q_pos as u32 + 1; // 1-indexed

        match pos {
            Some(0) => {
                // Ambiguously matching kmer
                // You can decide what to do here based on your specific requirements
            }
            Some(r_pos) => {
                let forward = r_pos > 0;
                let r_pos = r_pos.unsigned_abs(); // 1-indexed

                if let Some(a) = &mut alignment {
                    if a.forward == forward && a.r_end + 1 == r_pos && a.q_end + 1 == q_pos {
                        // Extend the current alignment
                        a.r_end = r_pos;
                        a.q_end = q_pos;
                    } else {
                        // Save the current alignment and start a new one
                        alignments.push(a.clone());
                        alignment = Some(Alignment {
                            q_start: q_pos,
                            q_end: q_pos,
                            r_start: r_pos,
                            r_end: r_pos,
                            forward,
                        });
                    }
                } else {
                    // Start a new alignment
                    alignment = Some(Alignment {
                        q_start: q_pos,
                        q_end: q_pos,
                        r_start: r_pos,
                        r_end: r_pos,
                        forward,
                    });
                }
            }
            None => {
                // Non-matching kmer
                if let Some(a) = alignment.take() {
                    // Save the current alignment if it exists
                    alignments.push(a);
                }
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
        let alignments = merge_segments(kmer_map);
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
        let alignments = merge_segments(kmer_map);
        assert_eq!(alignments.len(), 0);
    }

    #[test]
    fn test_ambiguous_kmers() {
        let kmer_map = vec![Some(0), Some(0), Some(0)];
        let alignments = merge_segments(kmer_map);
        // Assuming you don't create alignments for ambiguous kmers
        assert_eq!(alignments.len(), 0);
    }

    #[test]
    fn test_multiple_alignments() {
        let kmer_map = vec![Some(1), Some(2), None, Some(4), Some(5)];
        let alignments = merge_segments(kmer_map);
        assert_eq!(alignments.len(), 2);
        let a = &alignments[0];
        assert_eq!(a.q_start, 1);
        assert_eq!(a.q_end, 2);
        assert_eq!(a.r_start, 1);
        assert_eq!(a.r_end, 2);
        assert!(a.forward);
    }
}
