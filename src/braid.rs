//! Braid group representation.
//!
//! Implements B_n generators σ_i and the braid group relations:
//! - σ_i σ_j = σ_j σ_i for |i-j| ≥ 2 (far commutativity)
//! - σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1} (Yang-Baxter / braid relation)

use serde::{Serialize, Deserialize};
use nalgebra::{Complex, DMatrix};

/// A single braid generator σ_i (positive) or σ_i^{-1} (negative).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct BraidGenerator {
    pub index: usize,
    pub positive: bool,
}

impl BraidGenerator {
    pub fn sigma(index: usize) -> Self {
        Self { index, positive: true }
    }

    pub fn sigma_inv(index: usize) -> Self {
        Self { index, positive: false }
    }

    pub fn inverse(&self) -> Self {
        Self { index: self.index, positive: !self.positive }
    }
}

impl std::fmt::Display for BraidGenerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.positive {
            write!(f, "σ{}", self.index)
        } else {
            write!(f, "σ{}⁻¹", self.index)
        }
    }
}

/// A braid word: product of braid generators in B_n.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BraidWord {
    pub n: usize,
    pub generators: Vec<BraidGenerator>,
}

impl BraidWord {
    pub fn identity(n: usize) -> Self {
        Self { n, generators: vec![] }
    }

    pub fn sigma(&mut self, index: usize) -> &mut Self {
        self.generators.push(BraidGenerator::sigma(index));
        self
    }

    pub fn sigma_inv(&mut self, index: usize) -> &mut Self {
        self.generators.push(BraidGenerator::sigma_inv(index));
        self
    }

    pub fn inverse(&self) -> Self {
        let mut inv = Self {
            n: self.n,
            generators: Vec::with_capacity(self.generators.len()),
        };
        for gen in self.generators.iter().rev() {
            inv.generators.push(gen.inverse());
        }
        inv
    }

    pub fn is_valid(&self) -> bool {
        self.generators.iter().all(|g| g.index < self.n - 1)
    }

    pub fn compose(&self, other: &BraidWord) -> BraidWord {
        assert_eq!(self.n, other.n, "Braid words must have same number of strands");
        let mut result = self.clone();
        result.generators.extend(other.generators.iter().cloned());
        result
    }

    pub fn exponent_sum(&self) -> i64 {
        self.generators.iter().map(|g| if g.positive { 1i64 } else { -1i64 }).sum()
    }

    pub fn to_permutation(&self) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..self.n).collect();
        for gen in &self.generators {
            perm.swap(gen.index, gen.index + 1);
        }
        perm
    }

    pub fn is_identity_trivial(&self) -> bool {
        if self.exponent_sum() != 0 { return false; }
        let perm = self.to_permutation();
        perm.iter().enumerate().all(|(i, &p)| i == p)
    }
}

/// Build the unreduced Burau generator matrix for σ_i.
///
/// The unreduced Burau representation for B_n acts on C^n.
/// σ_i (0-indexed) is identity except at rows/columns i, i+1:
///   (i,i) = 1-t, (i,i+1) = 1, (i+1,i) = t, (i+1,i+1) = 0
fn burau_generator_matrix(n: usize, i: usize, t: Complex<f64>) -> DMatrix<Complex<f64>> {
    let one = Complex::new(1.0, 0.0);
    let one_minus_t = Complex::new(1.0, 0.0) - t;
    let mut mat = DMatrix::identity(n, n);

    mat[(i, i)] = one_minus_t;
    mat[(i, i + 1)] = one;
    mat[(i + 1, i)] = t;
    mat[(i + 1, i + 1)] = Complex::new(0.0, 0.0);

    mat
}

/// A braid with its Burau matrix representation.
#[derive(Debug, Clone)]
pub struct Braid {
    pub word: BraidWord,
}

impl Braid {
    pub fn new(word: BraidWord) -> Self {
        assert!(word.is_valid(), "Invalid braid word");
        Self { word }
    }

    pub fn identity(n: usize) -> Self {
        Self { word: BraidWord::identity(n) }
    }

    /// Compute the unreduced Burau representation matrix (n × n).
    pub fn burau_matrix(&self, t: Complex<f64>) -> DMatrix<Complex<f64>> {
        let n = self.word.n;
        if n <= 1 {
            return DMatrix::identity(1, 1);
        }
        let mut mat = DMatrix::identity(n, n);

        for gen in &self.word.generators {
            let gen_mat = burau_generator_matrix(n, gen.index, t);
            let inv = if !gen.positive {
                gen_mat.clone().try_inverse().unwrap_or(gen_mat.clone())
            } else {
                gen_mat
            };
            mat = inv * mat;
        }

        mat
    }

    /// Compute the reduced Burau representation matrix ((n-1) × (n-1)).
    /// Uses direct computation of the reduced representation.
    pub fn reduced_burau_matrix(&self, t: Complex<f64>) -> DMatrix<Complex<f64>> {
        let n = self.word.n;
        if n <= 1 {
            return DMatrix::identity(0, 0);
        }
        if n == 2 {
            let mut val = Complex::new(1.0, 0.0);
            for gen in &self.word.generators {
                let g = if gen.positive { -t } else { -Complex::new(1.0, 0.0) / t };
                val = g * val;
            }
            let mut result = DMatrix::zeros(1, 1);
            result[(0, 0)] = val;
            return result;
        }

        // For n >= 3: use the direct reduced Burau representation.
        // The reduced generators are (n-1)×(n-1) matrices:
        // σ_0: first 2×2 block is [[-t, 0], [-1, 1]], rest identity
        // σ_{n-2}: last 2×2 block is [[1, -t], [0, -t]], rest identity  
        // σ_i (0 < i < n-2): 3×3 block at (i-1, i, i+1) is [[1, 0, 0], [t, -t, 1], [0, 0, 1]], rest identity
        let k = n - 1;
        let mut mat = DMatrix::identity(k, k);

        for gen in &self.word.generators {
            let i = gen.index;
            let mut g = DMatrix::identity(k, k);

            if i == 0 {
                // σ_0 in reduced Burau
                g[(0, 0)] = -t;
                g[(0, 1)] = Complex::new(0.0, 0.0);
                g[(1, 0)] = Complex::new(-1.0, 0.0);
                g[(1, 1)] = Complex::new(1.0, 0.0);
            } else if i == k - 1 {
                // σ_{n-2} in reduced Burau
                g[(k-2, k-2)] = Complex::new(1.0, 0.0);
                g[(k-2, k-1)] = -t;
                g[(k-1, k-2)] = Complex::new(0.0, 0.0);
                g[(k-1, k-1)] = -t;
            } else {
                // Interior: 0 < i < k-1
                g[(i-1, i-1)] = Complex::new(1.0, 0.0);
                g[(i-1, i)] = Complex::new(0.0, 0.0);
                g[(i, i-1)] = t;
                g[(i, i)] = -t;
                g[(i, i+1)] = Complex::new(1.0, 0.0);
                g[(i+1, i)] = Complex::new(0.0, 0.0);
                g[(i+1, i+1)] = Complex::new(1.0, 0.0);
            }

            if gen.positive {
                mat = g * mat;
            } else {
                // Invert the generator matrix
                if let Some(g_inv) = g.clone().try_inverse() {
                    mat = g_inv * mat;
                } else {
                    mat = g * mat;
                }
            }
        }

        mat
    }

    /// Verify far commutativity: σ_i σ_j = σ_j σ_i for |i-j| ≥ 2.
    pub fn verify_far_commutativity(n: usize, i: usize, j: usize, t: Complex<f64>) -> bool {
        if i >= n - 1 || j >= n - 1 || (i as i64 - j as i64).unsigned_abs() < 2 {
            return true;
        }
        let mut bw1 = BraidWord::identity(n);
        bw1.sigma(i).sigma(j);
        let mut bw2 = BraidWord::identity(n);
        bw2.sigma(j).sigma(i);
        let b1 = Braid::new(bw1);
        let b2 = Braid::new(bw2);
        let m1 = b1.burau_matrix(t);
        let m2 = b2.burau_matrix(t);
        (m1 - m2).norm() < 1e-8
    }

    /// Verify Yang-Baxter: σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}.
    pub fn verify_yang_baxter(n: usize, i: usize, t: Complex<f64>) -> bool {
        if i >= n - 2 { return true; }
        let mut bw1 = BraidWord::identity(n);
        bw1.sigma(i).sigma(i + 1).sigma(i);
        let mut bw2 = BraidWord::identity(n);
        bw2.sigma(i + 1).sigma(i).sigma(i + 1);
        let b1 = Braid::new(bw1);
        let b2 = Braid::new(bw2);
        let m1 = b1.burau_matrix(t);
        let m2 = b2.burau_matrix(t);
        (m1 - m2).norm() < 1e-8
    }

    pub fn linking_number(&self) -> f64 {
        self.word.exponent_sum() as f64 / 2.0
    }
}

impl std::fmt::Display for Braid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "B{}(", self.word.n)?;
        for (i, gen) in self.word.generators.iter().enumerate() {
            if i > 0 { write!(f, " ")?; }
            write!(f, "{}", gen)?;
        }
        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_braid_identity() {
        let bw = BraidWord::identity(3);
        assert_eq!(bw.generators.len(), 0);
        assert!(bw.is_valid());
    }

    #[test]
    fn test_braid_generator_inverse() {
        let g = BraidGenerator::sigma(2);
        let gi = g.inverse();
        assert_eq!(gi.index, 2);
        assert!(!gi.positive);
        assert_eq!(gi.inverse(), g);
    }

    #[test]
    fn test_braid_word_inverse() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma(1).sigma_inv(0);
        let inv = bw.inverse();
        assert_eq!(inv.generators.len(), 3);
        assert_eq!(inv.generators[0], BraidGenerator::sigma(0));
        assert_eq!(inv.generators[1], BraidGenerator::sigma_inv(1));
        assert_eq!(inv.generators[2], BraidGenerator::sigma_inv(0));
    }

    #[test]
    fn test_braid_validity() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma(1);
        assert!(bw.is_valid());

        let bw2 = BraidWord { n: 3, generators: vec![BraidGenerator::sigma(2)] };
        assert!(!bw2.is_valid());
    }

    #[test]
    fn test_exponent_sum() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma(1).sigma_inv(0);
        assert_eq!(bw.exponent_sum(), 1);
    }

    #[test]
    fn test_permutation() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0);
        let perm = bw.to_permutation();
        assert_eq!(perm, vec![1, 0, 2]);
    }

    #[test]
    fn test_compose() {
        let mut bw1 = BraidWord::identity(3);
        bw1.sigma(0);
        let mut bw2 = BraidWord::identity(3);
        bw2.sigma(1);
        let composed = bw1.compose(&bw2);
        assert_eq!(composed.generators.len(), 2);
    }

    #[test]
    fn test_yang_baxter_relation() {
        let t = Complex::new(0.5, 0.3);
        assert!(Braid::verify_yang_baxter(4, 0, t));
        assert!(Braid::verify_yang_baxter(4, 1, t));
        assert!(Braid::verify_yang_baxter(5, 2, t));
    }

    #[test]
    fn test_far_commutativity() {
        let t = Complex::new(0.7, 0.2);
        assert!(Braid::verify_far_commutativity(5, 0, 2, t));
        assert!(Braid::verify_far_commutativity(5, 0, 3, t));
        assert!(Braid::verify_far_commutativity(6, 1, 3, t));
    }

    #[test]
    fn test_burau_identity() {
        let t = Complex::new(0.6, 0.0);
        let b = Braid::identity(3);
        let mat = b.burau_matrix(t);
        let identity = DMatrix::identity(3, 3);
        assert!((mat - identity).norm() < 1e-10);
    }

    #[test]
    fn test_braid_display() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma(1);
        let b = Braid::new(bw);
        let s = format!("{}", b);
        assert!(s.contains("σ0"));
        assert!(s.contains("σ1"));
    }

    #[test]
    fn test_linking_number() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma(0);
        let b = Braid::new(bw);
        assert_eq!(b.linking_number(), 1.0);
    }

    #[test]
    fn test_identity_trivial() {
        let mut bw = BraidWord::identity(3);
        bw.sigma(0).sigma_inv(0);
        assert!(bw.is_identity_trivial());

        let mut bw2 = BraidWord::identity(3);
        bw2.sigma(0);
        assert!(!bw2.is_identity_trivial());
    }
}
