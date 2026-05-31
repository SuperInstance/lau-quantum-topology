//! Modular Tensor Categories (MTC).
//!
//! Implements the algebraic structure of modular tensor categories:
//! - Simple objects and fusion rules
//! - S-matrix and T-matrix (modular data)
//! - Quantum dimensions from fusion rules
//! - Verlinde formula for fusion coefficients

use serde::{Serialize, Deserialize};
use nalgebra::{Complex, DMatrix};

/// A simple object in a modular tensor category.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SimpleObject {
    /// Label/name of the simple object.
    pub label: String,
    /// Quantum dimension.
    pub quantum_dimension: f64,
    /// Topological spin (as a rational fraction p/q of 2π).
    pub twist_fraction: (i64, i64),
    /// Whether this is the unit object.
    pub is_unit: bool,
}

impl SimpleObject {
    /// Create a new simple object.
    pub fn new(label: impl Into<String>, quantum_dimension: f64, twist_fraction: (i64, i64)) -> Self {
        Self {
            label: label.into(),
            quantum_dimension,
            twist_fraction,
            is_unit: false,
        }
    }

    /// The unit/vacuum object.
    pub fn unit() -> Self {
        Self {
            label: "1".to_string(),
            quantum_dimension: 1.0,
            twist_fraction: (0, 1),
            is_unit: true,
        }
    }

    /// Compute the topological twist θ = e^{2πi * p/q}.
    pub fn twist(&self) -> Complex<f64> {
        let (p, q) = self.twist_fraction;
        let angle = 2.0 * std::f64::consts::PI * (p as f64) / (q as f64);
        Complex::new(angle.cos(), angle.sin())
    }
}

/// A fusion rule: a × b = Σ N_{ab}^c c.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FusionRule {
    pub a: String,
    pub b: String,
    /// (c, multiplicity) pairs.
    pub outcomes: Vec<(String, usize)>,
}

impl FusionRule {
    /// Create a simple fusion rule a × b → c with multiplicity 1.
    pub fn simple(a: impl Into<String>, b: impl Into<String>, c: impl Into<String>) -> Self {
        Self {
            a: a.into(),
            b: b.into(),
            outcomes: vec![(c.into(), 1)],
        }
    }

    /// Create a fusion rule with multiple outcomes.
    pub fn multi(a: impl Into<String>, b: impl Into<String>, outcomes: Vec<(String, usize)>) -> Self {
        Self { a: a.into(), b: b.into(), outcomes }
    }
}

/// A modular tensor category.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModularTensorCategory {
    /// Name of the MTC.
    pub name: String,
    /// Simple objects.
    pub simple_objects: Vec<SimpleObject>,
    /// Fusion rules.
    pub fusion_rules: Vec<FusionRule>,
    /// S-matrix (modular transformation).
    pub s_matrix: Vec<Vec<Complex<f64>>>,
    /// T-matrix (twist).
    pub t_matrix: Vec<Vec<Complex<f64>>>,
}

impl ModularTensorCategory {
    /// Construct the Fibonacci (Fib) modular tensor category.
    /// Simple objects: {1, τ}
    /// Fusion: τ × τ = 1 + τ
    pub fn fibonacci() -> Self {
        let simple_objects = vec![
            SimpleObject::unit(),
            SimpleObject::new("τ", (1.0 + 5.0_f64.sqrt()) / 2.0, (2, 5)),
        ];

        let fusion_rules = vec![
            FusionRule::simple("1", "1", "1"),
            FusionRule::simple("1", "τ", "τ"),
            FusionRule::simple("τ", "1", "τ"),
            FusionRule::multi("τ", "τ", vec![
                ("1".to_string(), 1),
                ("τ".to_string(), 1),
            ]),
        ];

        // S-matrix for Fibonacci MTC (unitary normalization):
        // S = (1/D) * [[1, φ], [φ, -1]]
        // where D² = 1 + φ² and φ = (1+√5)/2
        // This satisfies S†S = I and Verlinde formula.
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let d_sq = 1.0 + phi * phi;
        let d = d_sq.sqrt();
        let s11 = -1.0 / d;  // S_{τ,τ}
        let s_matrix = vec![
            vec![Complex::new(1.0 / d, 0.0), Complex::new(phi / d, 0.0)],
            vec![Complex::new(phi / d, 0.0), Complex::new(s11, 0.0)],
        ];

        // T-matrix: diagonal with twists
        // θ_1 = 1, θ_τ = e^{4πi/5}
        // But for the Fibonacci MTC, the correct twist is θ_τ = e^{4πi/5}
        // However, the modular relation (ST)^3 = S² requires specific twist values.
        // For Fibonacci: θ_τ = e^{4πi/5}
        let tau_angle = 4.0 * std::f64::consts::PI / 5.0;
        let t_matrix = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            vec![Complex::new(0.0, 0.0), Complex::new(tau_angle.cos(), tau_angle.sin())],
        ];

        Self {
            name: "Fibonacci (Fib)".to_string(),
            simple_objects,
            fusion_rules,
            s_matrix,
            t_matrix,
        }
    }

    /// Construct the Ising modular tensor category.
    /// Simple objects: {1, ψ, σ}
    /// Fusion: σ × σ = 1 + ψ, σ × ψ = σ, ψ × ψ = 1
    pub fn ising() -> Self {
        let simple_objects = vec![
            SimpleObject::unit(),
            SimpleObject::new("ψ", 1.0, (1, 2)),
            SimpleObject::new("σ", 2.0_f64.sqrt(), (1, 16)),
        ];

        let fusion_rules = vec![
            FusionRule::simple("1", "1", "1"),
            FusionRule::simple("1", "ψ", "ψ"),
            FusionRule::simple("ψ", "1", "ψ"),
            FusionRule::simple("1", "σ", "σ"),
            FusionRule::simple("σ", "1", "σ"),
            FusionRule::simple("ψ", "ψ", "1"),
            FusionRule::simple("σ", "ψ", "σ"),
            FusionRule::simple("ψ", "σ", "σ"),
            FusionRule::multi("σ", "σ", vec![
                ("1".to_string(), 1),
                ("ψ".to_string(), 1),
            ]),
        ];

        // S-matrix for Ising:
        // S = (1/2) * [[1, 1, √2], [1, 1, -√2], [√2, -√2, 0]]
        let inv2 = 0.5;
        let sqrt2 = 2.0_f64.sqrt();
        let s_matrix = vec![
            vec![Complex::new(inv2, 0.0), Complex::new(inv2, 0.0), Complex::new(sqrt2 * inv2, 0.0)],
            vec![Complex::new(inv2, 0.0), Complex::new(inv2, 0.0), Complex::new(-sqrt2 * inv2, 0.0)],
            vec![Complex::new(sqrt2 * inv2, 0.0), Complex::new(-sqrt2 * inv2, 0.0), Complex::new(0.0, 0.0)],
        ];

        // T-matrix: θ_1 = 1, θ_ψ = -1, θ_σ = e^{iπ/8}
        let sigma_angle = std::f64::consts::PI / 8.0;
        let t_matrix = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            vec![Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0), Complex::new(0.0, 0.0)],
            vec![Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(sigma_angle.cos(), sigma_angle.sin())],
        ];

        Self {
            name: "Ising".to_string(),
            simple_objects,
            fusion_rules,
            s_matrix,
            t_matrix,
        }
    }

    /// Get the total quantum dimension D = √(Σ d_i²).
    pub fn total_quantum_dimension(&self) -> f64 {
        let sum: f64 = self.simple_objects.iter().map(|o| o.quantum_dimension.powi(2)).sum();
        sum.sqrt()
    }

    /// Get the S-matrix as a nalgebra matrix.
    pub fn s_matrix_na(&self) -> DMatrix<Complex<f64>> {
        let n = self.s_matrix.len();
        let mut vals = Vec::with_capacity(n * n);
        for row in &self.s_matrix {
            for val in row {
                vals.push(*val);
            }
        }
        DMatrix::from_row_slice(n, n, &vals)
    }

    /// Get the T-matrix as a nalgebra matrix.
    pub fn t_matrix_na(&self) -> DMatrix<Complex<f64>> {
        let n = self.t_matrix.len();
        let mut vals = Vec::with_capacity(n * n);
        for row in &self.t_matrix {
            for val in row {
                vals.push(*val);
            }
        }
        DMatrix::from_row_slice(n, n, &vals)
    }

    /// Verify S-matrix unitarity: S† S = I.
    pub fn verify_s_unitarity(&self) -> bool {
        let s = self.s_matrix_na();
        let n = s.nrows();
        let product = s.adjoint() * &s;
        let identity = DMatrix::identity(n, n);
        (product - identity).norm() < 1e-8
    }

    /// Verify the Verlinde formula:
    /// N_{ij}^k = Σ_l (S_{il} S_{jl} S_{kl}^*) / S_{0l}
    pub fn verify_verlinde(&self) -> bool {
        let s = self.s_matrix_na();
        let n = s.nrows();

        let labels: Vec<&str> = self.simple_objects.iter().map(|o| o.label.as_str()).collect();

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    // Compute Verlinde formula: N_{ij}^k = Σ_l S_{il} S_{jl} conj(S_{kl}) / S_{0l}
                    let mut n_ijk = Complex::new(0.0, 0.0);
                    for l in 0..n {
                        let s_0l = s[(0, l)];
                        if s_0l.norm() < 1e-10 {
                            continue;
                        }
                        n_ijk += s[(i, l)] * s[(j, l)] * s[(k, l)].conj() / s_0l;
                    }

                    // Compare with actual fusion rules
                    let actual = self.get_fusion_coeff(labels[i], labels[j], labels[k]);
                    // The Verlinde formula gives the actual coefficient as a real number
                    // (should be a non-negative integer)
                    let diff = (n_ijk.re - actual as f64).abs() + n_ijk.im.abs();
                    if diff > 1e-4 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Get the fusion coefficient N_{ab}^c.
    fn get_fusion_coeff(&self, a: &str, b: &str, c: &str) -> usize {
        for rule in &self.fusion_rules {
            if rule.a == a && rule.b == b {
                for (outcome, mult) in &rule.outcomes {
                    if outcome == c {
                        return *mult;
                    }
                }
            }
        }
        0
    }

    /// Verify the modular relation: (ST)^3 = S².
    /// In a modular tensor category, (ST)^3 should equal a scalar multiple of S².
    pub fn verify_modular_relation(&self) -> bool {
        let s = self.s_matrix_na();
        let t = self.t_matrix_na();

        // (ST)^3
        let st = &s * &t;
        let st3 = &st * &st * &st;

        // S²  
        let s2 = &s * &s;

        // Check if (ST)^3 is proportional to S²
        // They should be equal up to a global phase for a valid MTC
        // More precisely: S^2 = C (charge conjugation) and (ST)^3 = S^2 * C
        // So (ST)^3 should equal some multiple of S^2
        if s2.norm() < 1e-10 {
            return true;
        }
        // Find the ratio of any non-zero element
        let mut ratio = None;
        for i in 0..s2.nrows() {
            for j in 0..s2.ncols() {
                if s2[(i, j)].norm() > 1e-8 {
                    ratio = Some(st3[(i, j)] / s2[(i, j)]);
                    break;
                }
            }
            if ratio.is_some() { break; }
        }

        match ratio {
            Some(r) => {
                // Check that st3 = r * s2 element-wise
                for i in 0..s2.nrows() {
                    for j in 0..s2.ncols() {
                        let expected = r * s2[(i, j)];
                        if (st3[(i, j)] - expected).norm() > 1e-4 {
                            return false;
                        }
                    }
                }
                // The ratio should have unit modulus
                (r.norm() - 1.0).abs() < 1e-4
            }
            None => true,
        }
    }

    /// Verify T-matrix: should be diagonal with unit modulus entries.
    pub fn verify_t_diagonal(&self) -> bool {
        let t = self.t_matrix_na();
        let n = t.nrows();
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    if (t[(i, i)].norm() - 1.0).abs() > 1e-8 {
                        return false;
                    }
                } else if t[(i, j)].norm() > 1e-8 {
                    return false;
                }
            }
        }
        true
    }

    /// Compute quantum dimensions from the S-matrix: d_i = S_{0i} / S_{00}.
    pub fn quantum_dimensions_from_s(&self) -> Vec<f64> {
        let s = self.s_matrix_na();
        let n = s.nrows();
        let s00 = s[(0, 0)];
        (0..n).map(|i| (s[(0, i)] / s00).re).collect()
    }

    /// Compute the Gauss sum: Σ_i d_i² θ_i
    pub fn gauss_sum(&self) -> Complex<f64> {
        self.simple_objects
            .iter()
            .map(|o| Complex::new(o.quantum_dimension.powi(2), 0.0) * o.twist())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fibonacci_s_unitarity() {
        let mtc = ModularTensorCategory::fibonacci();
        assert!(mtc.verify_s_unitarity());
    }

    #[test]
    fn test_ising_s_unitarity() {
        let mtc = ModularTensorCategory::ising();
        assert!(mtc.verify_s_unitarity());
    }

    #[test]
    fn test_fibonacci_total_quantum_dimension() {
        let mtc = ModularTensorCategory::fibonacci();
        let d = mtc.total_quantum_dimension();
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let expected = (1.0 + phi * phi).sqrt();
        assert!((d - expected).abs() < 1e-8);
    }

    #[test]
    fn test_ising_total_quantum_dimension() {
        let mtc = ModularTensorCategory::ising();
        let d = mtc.total_quantum_dimension();
        // d² = 1 + 1 + 2 = 4, so d = 2
        assert!((d - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_fibonacci_fusion() {
        let mtc = ModularTensorCategory::fibonacci();
        assert_eq!(mtc.get_fusion_coeff("τ", "τ", "1"), 1);
        assert_eq!(mtc.get_fusion_coeff("τ", "τ", "τ"), 1);
        assert_eq!(mtc.get_fusion_coeff("1", "τ", "τ"), 1);
    }

    #[test]
    fn test_ising_fusion() {
        let mtc = ModularTensorCategory::ising();
        assert_eq!(mtc.get_fusion_coeff("σ", "σ", "1"), 1);
        assert_eq!(mtc.get_fusion_coeff("σ", "σ", "ψ"), 1);
        assert_eq!(mtc.get_fusion_coeff("ψ", "ψ", "1"), 1);
        assert_eq!(mtc.get_fusion_coeff("σ", "ψ", "σ"), 1);
    }

    #[test]
    fn test_fibonacci_verlinde() {
        let mtc = ModularTensorCategory::fibonacci();
        assert!(mtc.verify_verlinde());
    }

    #[test]
    fn test_ising_verlinde() {
        let mtc = ModularTensorCategory::ising();
        assert!(mtc.verify_verlinde());
    }

    #[test]
    fn test_fibonacci_t_diagonal() {
        let mtc = ModularTensorCategory::fibonacci();
        assert!(mtc.verify_t_diagonal());
    }

    #[test]
    fn test_ising_t_diagonal() {
        let mtc = ModularTensorCategory::ising();
        assert!(mtc.verify_t_diagonal());
    }

    #[test]
    fn test_quantum_dimensions_from_s_fibonacci() {
        let mtc = ModularTensorCategory::fibonacci();
        let dims = mtc.quantum_dimensions_from_s();
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((dims[0] - 1.0).abs() < 1e-8);
        assert!((dims[1] - phi).abs() < 1e-6);
    }

    #[test]
    fn test_quantum_dimensions_from_s_ising() {
        let mtc = ModularTensorCategory::ising();
        let dims = mtc.quantum_dimensions_from_s();
        assert!((dims[0] - 1.0).abs() < 1e-6);
        assert!((dims[1] - 1.0).abs() < 1e-6);
        assert!((dims[2] - 2.0_f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_simple_object_twist() {
        let unit = SimpleObject::unit();
        let twist = unit.twist();
        assert!((twist.re - 1.0).abs() < 1e-10);
        assert!(twist.im.abs() < 1e-10);
    }

    #[test]
    fn test_fibonacci_modular_relation() {
        let mtc = ModularTensorCategory::fibonacci();
        assert!(mtc.verify_modular_relation());
    }

    #[test]
    fn test_ising_modular_relation() {
        let mtc = ModularTensorCategory::ising();
        assert!(mtc.verify_modular_relation());
    }

    #[test]
    fn test_gauss_sum_fibonacci() {
        let mtc = ModularTensorCategory::fibonacci();
        let gs = mtc.gauss_sum();
        // Gauss sum should have modulus equal to D
        let d = mtc.total_quantum_dimension();
        assert!((gs.norm() - d).abs() < 1e-6);
    }

    #[test]
    fn test_gauss_sum_ising() {
        let mtc = ModularTensorCategory::ising();
        let gs = mtc.gauss_sum();
        let d = mtc.total_quantum_dimension();
        assert!((gs.norm() - d).abs() < 1e-6);
    }
}

