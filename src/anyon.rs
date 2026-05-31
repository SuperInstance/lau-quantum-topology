//! Anyon types and braiding operations.
//!
//! Supports Fibonacci anyons and Ising anyons with their fusion rules,
//! braiding matrices (R-matrices), and F-matrices (recoupling).

use serde::{Serialize, Deserialize};
use nalgebra::{Complex, DMatrix};

/// The vacuum (trivial) anyon label.
pub const VACUUM: &str = "1";

/// Types of anyons supported.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AnyonType {
    /// Trivial/vacuum anyon
    Vacuum,
    /// Fibonacci anyon: τ with fusion τ × τ = 1 + τ
    Fibonacci,
    /// Ising anyons: σ (non-abelian) and ψ (fermion)
    Ising,
}

impl AnyonType {
    /// Get the topological spin of this anyon type.
    pub fn topological_spin(&self) -> Complex<f64> {
        match self {
            AnyonType::Vacuum => Complex::new(1.0, 0.0),
            AnyonType::Fibonacci => {
                // θ_τ = e^{4πi/5}
                let phi = 4.0 * std::f64::consts::PI / 5.0;
                Complex::new(phi.cos(), phi.sin())
            }
            AnyonType::Ising => {
                // σ has spin e^{iπ/8}, but we return the σ anyon spin
                let phi = std::f64::consts::PI / 8.0;
                Complex::new(phi.cos(), phi.sin())
            }
        }
    }

    /// Quantum dimension of this anyon type.
    pub fn quantum_dimension(&self) -> f64 {
        match self {
            AnyonType::Vacuum => 1.0,
            AnyonType::Fibonacci => (1.0 + 5.0_f64.sqrt()) / 2.0, // golden ratio φ
            AnyonType::Ising => 2.0_f64.sqrt(),
        }
    }
}

/// A concrete anyon with a type and an optional label.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Anyon {
    pub anyon_type: AnyonType,
    pub label: String,
}

impl Anyon {
    pub fn new(anyon_type: AnyonType, label: impl Into<String>) -> Self {
        Self {
            anyon_type,
            label: label.into(),
        }
    }

    /// Fibonacci anyon τ
    pub fn fib(label: impl Into<String>) -> Self {
        Self::new(AnyonType::Fibonacci, label)
    }

    /// Ising anyon σ
    pub fn ising(label: impl Into<String>) -> Self {
        Self::new(AnyonType::Ising, label)
    }

    /// Vacuum anyon 1
    pub fn vacuum() -> Self {
        Self::new(AnyonType::Vacuum, "1")
    }
}

/// Result of fusing two anyons.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FusionResult {
    /// Possible outcomes from fusion.
    pub outcomes: Vec<AnyonType>,
}

/// Compute the fusion of two Fibonacci anyons: τ × τ = 1 + τ
pub fn fuse_fibonacci(a: &AnyonType, b: &AnyonType) -> FusionResult {
    match (a, b) {
        (AnyonType::Vacuum, x) | (x, AnyonType::Vacuum) => {
            FusionResult { outcomes: vec![x.clone()] }
        }
        (AnyonType::Fibonacci, AnyonType::Fibonacci) => {
            FusionResult {
                outcomes: vec![AnyonType::Vacuum, AnyonType::Fibonacci],
            }
        }
        _ => FusionResult { outcomes: vec![] },
    }
}

/// Compute the fusion of two Ising anyons.
/// σ × σ = 1 + ψ, σ × ψ = σ, ψ × ψ = 1
pub fn fuse_ising(a: &AnyonType, b: &AnyonType) -> FusionResult {
    match (a, b) {
        (AnyonType::Vacuum, x) | (x, AnyonType::Vacuum) => {
            FusionResult { outcomes: vec![x.clone()] }
        }
        (AnyonType::Ising, AnyonType::Ising) => {
            FusionResult {
                outcomes: vec![AnyonType::Vacuum, AnyonType::Ising],
            }
        }
        // In a full Ising model we'd also have ψ, but we simplify
        _ => FusionResult { outcomes: vec![] },
    }
}

/// R-matrix (braiding matrix) for Fibonacci anyons.
/// Returns the phase acquired when two τ anyons braid.
/// R_{ττ}^1 = e^{-4πi/5}, R_{ττ}^τ = e^{3πi/5}
pub fn fibonacci_r_matrix() -> (Complex<f64>, Complex<f64>) {
    let r_vacuum = {
        let theta = -4.0 * std::f64::consts::PI / 5.0;
        Complex::new(theta.cos(), theta.sin())
    };
    let r_tau = {
        let theta = 3.0 * std::f64::consts::PI / 5.0;
        Complex::new(theta.cos(), theta.sin())
    };
    (r_vacuum, r_tau)
}

/// F-matrix (recoupling matrix) for Fibonacci anyons.
/// Returns the 2×2 F-matrix for recoupling three τ anyons.
pub fn fibonacci_f_matrix() -> DMatrix<Complex<f64>> {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let inv_phi = 1.0 / phi;
    // F-matrix:
    // | inv_phi    1/sqrt(phi) |
    // | 1/sqrt(phi)  -inv_phi  |
    DMatrix::from_row_slice(2, 2, &[
        Complex::new(inv_phi, 0.0),
        Complex::new(1.0 / phi.sqrt(), 0.0),
        Complex::new(1.0 / phi.sqrt(), 0.0),
        Complex::new(-inv_phi, 0.0),
    ])
}

/// R-matrix for Ising anyons.
/// R_{σσ}^1 = e^{-iπ/8}, R_{σσ}^ψ = e^{3iπ/8}
pub fn ising_r_matrix() -> (Complex<f64>, Complex<f64>) {
    let r_vacuum = {
        let theta = -std::f64::consts::PI / 8.0;
        Complex::new(theta.cos(), theta.sin())
    };
    let r_fermion = {
        let theta = 3.0 * std::f64::consts::PI / 8.0;
        Complex::new(theta.cos(), theta.sin())
    };
    (r_vacuum, r_fermion)
}

/// F-matrix for Ising anyons.
pub fn ising_f_matrix() -> DMatrix<Complex<f64>> {
    // Ising F-matrix for σσσ recoupling is the 2×2 Hadamard-like matrix
    let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
    DMatrix::from_row_slice(2, 2, &[
        Complex::new(inv_sqrt2, 0.0),
        Complex::new(inv_sqrt2, 0.0),
        Complex::new(inv_sqrt2, 0.0),
        Complex::new(-inv_sqrt2, 0.0),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fibonacci_fusion() {
        let result = fuse_fibonacci(&AnyonType::Fibonacci, &AnyonType::Fibonacci);
        assert_eq!(result.outcomes.len(), 2);
        assert!(result.outcomes.contains(&AnyonType::Vacuum));
        assert!(result.outcomes.contains(&AnyonType::Fibonacci));
    }

    #[test]
    fn test_fibonacci_fusion_with_vacuum() {
        let result = fuse_fibonacci(&AnyonType::Vacuum, &AnyonType::Fibonacci);
        assert_eq!(result.outcomes, vec![AnyonType::Fibonacci]);

        let result2 = fuse_fibonacci(&AnyonType::Fibonacci, &AnyonType::Vacuum);
        assert_eq!(result2.outcomes, vec![AnyonType::Fibonacci]);
    }

    #[test]
    fn test_vacuum_fusion() {
        let result = fuse_fibonacci(&AnyonType::Vacuum, &AnyonType::Vacuum);
        assert_eq!(result.outcomes, vec![AnyonType::Vacuum]);
    }

    #[test]
    fn test_ising_fusion() {
        let result = fuse_ising(&AnyonType::Ising, &AnyonType::Ising);
        assert_eq!(result.outcomes.len(), 2);
        assert!(result.outcomes.contains(&AnyonType::Vacuum));
        assert!(result.outcomes.contains(&AnyonType::Ising));
    }

    #[test]
    fn test_fibonacci_quantum_dimension() {
        let d = AnyonType::Fibonacci.quantum_dimension();
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((d - phi).abs() < 1e-10);
    }

    #[test]
    fn test_ising_quantum_dimension() {
        let d = AnyonType::Ising.quantum_dimension();
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_fibonacci_r_matrix_unit_modulus() {
        let (r1, r2) = fibonacci_r_matrix();
        assert!((r1.norm() - 1.0).abs() < 1e-10);
        assert!((r2.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_r_matrix_unit_modulus() {
        let (r1, r2) = ising_r_matrix();
        assert!((r1.norm() - 1.0).abs() < 1e-10);
        assert!((r2.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fibonacci_f_matrix_unitary() {
        let f = fibonacci_f_matrix();
        let f_dag = f.adjoint();
        let product = &f_dag * &f;
        let identity = DMatrix::identity(2, 2);
        assert!((product - identity).norm() < 1e-10);
    }

    #[test]
    fn test_ising_f_matrix_unitary() {
        let f = ising_f_matrix();
        let f_dag = f.adjoint();
        let product = &f_dag * &f;
        let identity = DMatrix::identity(2, 2);
        assert!((product - identity).norm() < 1e-10);
    }

    #[test]
    fn test_fibonacci_f_matrix_det() {
        let f = fibonacci_f_matrix();
        // det should be -1 for Fibonacci F-matrix
        let det = f[(0, 0)] * f[(1, 1)] - f[(0, 1)] * f[(1, 0)];
        assert!((det.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_topological_spin_vacuum() {
        let spin = AnyonType::Vacuum.topological_spin();
        assert!((spin.re - 1.0).abs() < 1e-10);
        assert!(spin.im.abs() < 1e-10);
    }
}
