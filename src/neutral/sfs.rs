use crate::model::{CountRow, SFS_PROFILE};
use rsomics_common::{Result, RsomicsError};
use serde_json::{Value, json};

pub(crate) struct Spectrum {
    names: [String; 2],
    chromosomes: [usize; 2],
    folded: [Vec<u64>; 2],
    unfolded: [Vec<u64>; 2],
    joint_folded: Vec<Vec<u64>>,
    joint_unfolded: Vec<Vec<u64>>,
    eligibility: [[u64; 6]; 2],
    joint: [u64; 4],
}

impl Spectrum {
    pub fn new(names: [String; 2], sizes: [usize; 2], budget: usize) -> Result<Self> {
        let chromosomes = sizes.map(|value| 2 * value);
        let cells = (chromosomes[0] + 1)
            .checked_mul(chromosomes[1] + 1)
            .ok_or_else(|| RsomicsError::InvalidInput("SFS axes overflow".into()))?;
        if cells.checked_mul(16).is_none_or(|value| value > budget / 4) {
            return Err(RsomicsError::InvalidInput(
                "dense SFS axes exceed the planning budget".into(),
            ));
        }
        Ok(Self {
            names,
            chromosomes,
            folded: chromosomes.map(|n| vec![0; n / 2 + 1]),
            unfolded: chromosomes.map(|n| vec![0; n + 1]),
            joint_folded: vec![vec![0; chromosomes[1] / 2 + 1]; chromosomes[0] / 2 + 1],
            joint_unfolded: vec![vec![0; chromosomes[1] + 1]; chromosomes[0] + 1],
            eligibility: [[0; 6]; 2],
            joint: [0; 4],
        })
    }

    pub fn push(&mut self, row: &CountRow, ancestor: i32) {
        let biallelic = row.alleles.len() == 2;
        let complete = std::array::from_fn::<_, 2, _>(|pop| {
            biallelic
                && row.populations[pop]
                    .copies
                    .iter()
                    .map(|&value| value as usize)
                    .sum::<usize>()
                    == self.chromosomes[pop]
        });
        let usable = (0..2).contains(&ancestor);
        for (pop, &is_complete) in complete.iter().enumerate() {
            let eligibility = &mut self.eligibility[pop];
            eligibility[0] += 1;
            eligibility[1] += u64::from(biallelic);
            if is_complete {
                let ac = &row.populations[pop].copies;
                self.folded[pop][ac[0].min(ac[1]) as usize] += 1;
                eligibility[2] += 1;
                if usable {
                    self.unfolded[pop][ac[1 - ancestor as usize] as usize] += 1;
                    eligibility[3] += 1;
                } else {
                    eligibility[5] += 1;
                }
            } else if biallelic {
                eligibility[4] += 1;
            }
        }
        if complete.iter().all(|&value| value) {
            let first = &row.populations[0].copies;
            let second = &row.populations[1].copies;
            self.joint_folded[first[0].min(first[1]) as usize]
                [second[0].min(second[1]) as usize] += 1;
            self.joint[0] += 1;
            if usable {
                self.joint_unfolded[first[1 - ancestor as usize] as usize]
                    [second[1 - ancestor as usize] as usize] += 1;
                self.joint[1] += 1;
            } else {
                self.joint[3] += 1;
            }
        } else if biallelic {
            self.joint[2] += 1;
        }
    }

    pub fn finish(self) -> Result<Value> {
        if self
            .eligibility
            .iter()
            .any(|value| value[2] == 0 || value[3] == 0)
            || self.joint[0] == 0
            || self.joint[1] == 0
        {
            return Err(RsomicsError::InvalidInput(
                "complete workflow requires nonempty folded and polarized SFS eligibility".into(),
            ));
        }
        let mut populations = serde_json::Map::new();
        for pop in 0..2 {
            let e = self.eligibility[pop];
            populations.insert(self.names[pop].clone(), json!({"folded":self.folded[pop],"unfolded":self.unfolded[pop],
                "eligibility":{"retained":e[0],"biallelic":e[1],"folded":e[2],"unfolded":e[3],"incomplete_biallelic":e[4],"unusable_ancestral_complete":e[5]}}));
        }
        let e = self.joint;
        Ok(
            json!({"profile":SFS_PROFILE,"chromosome_counts":{&self.names[0]:self.chromosomes[0],&self.names[1]:self.chromosomes[1]},"populations":populations,
            "joint":{"populations":self.names,"folded":self.joint_folded,"unfolded":self.joint_unfolded,
                "eligibility":{"folded":e[0],"unfolded":e[1],"incomplete_biallelic":e[2],"unusable_ancestral_complete":e[3]}}}),
        )
    }
}
