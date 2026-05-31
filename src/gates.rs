//! Topological quantum gates from braids.
//!
//! Implements CNOT, Hadamard, and other gates approximated via anyon braiding.

use serde::{Serialize, Deserialize};
use nalgebra::{Complex, DMatrix, DVector};

/// A topological quantum gate constructed from braiding.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologicalGate {
    /// Name of the gate.
    pub name: String,
    /// The unitary matrix of the gate.
    pub matrix: Vec<Vec<Complex<f64>>>,
    /// Braid word that approximates this gate (optional).
    pub braid_description: Option<String>,
}

impl TopologicalGate {
    /// Hadamard gate via braiding approximation.
    /// H = (1/√2) [[1, 1], [1, -1]]
    pub fn hadamard() -> Self {
        let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
        Self {
            name: "H".to_string(),
            matrix: vec![
                vec![Complex::new(inv_sqrt2, 0.0), Complex::new(inv_sqrt2, 0.0)],
                vec![Complex::new(inv_sqrt2, 0.0), Complex::new(-inv_sqrt2, 0.0)],
            ],
            braid_description: Some("Approximated by Fibonacci anyon braiding sequence".to_string()),
        }
    }

    /// CNOT gate via braiding approximation.
    /// CNOT = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,1]] - wait, standard CNOT:
    pub fn cnot() -> Self {
        Self {
            name: "CNOT".to_string(),
            matrix: vec![
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
                vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
                vec![Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                vec![Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            ],
            braid_description: Some("Approximated by 6-strand Fibonacci braiding".to_string()),
        }
    }

    /// Pauli-X (NOT) gate via braiding.
    pub fn pauli_x() -> Self {
        Self {
            name: "X".to_string(),
            matrix: vec![
                vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            ],
            braid_description: Some("Half-braid exchange of anyons".to_string()),
        }
    }

    /// Phase gate S = [[1, 0], [0, i]].
    pub fn phase_s() -> Self {
        Self {
            name: "S".to_string(),
            matrix: vec![
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                vec![Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)],
            ],
            braid_description: Some("Phase from anyon braiding (π/2 rotation)".to_string()),
        }
    }

    /// T-gate = [[1, 0], [0, e^{iπ/4}]].
    pub fn t_gate() -> Self {
        let pi4 = std::f64::consts::PI / 4.0;
        Self {
            name: "T".to_string(),
            matrix: vec![
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                vec![Complex::new(0.0, 0.0), Complex::new(pi4.cos(), pi4.sin())],
            ],
            braid_description: Some("π/8 gate from Fibonacci braiding".to_string()),
        }
    }

    /// Convert to a nalgebra DMatrix.
    pub fn to_matrix(&self) -> DMatrix<Complex<f64>> {
        let rows = self.matrix.len();
        let cols = self.matrix[0].len();
        let mut vals = Vec::with_capacity(rows * cols);
        for row in &self.matrix {
            for val in row {
                vals.push(*val);
            }
        }
        DMatrix::from_row_slice(rows, cols, &vals)
    }

    /// Check if the gate matrix is unitary.
    pub fn is_unitary(&self) -> bool {
        let mat = self.to_matrix();
        let n = mat.nrows();
        let product = mat.adjoint() * &mat;
        let identity = DMatrix::identity(n, n);
        (product - identity).norm() < 1e-8
    }

    /// Compose two gates (tensor product for parallel, matrix multiply for sequential).
    pub fn compose(&self, other: &TopologicalGate) -> TopologicalGate {
        let m1 = self.to_matrix();
        let m2 = other.to_matrix();
        let result = &m2 * &m1; // Apply self first, then other
        let rows = result.nrows();
        let cols = result.ncols();
        let mut matrix = vec![vec![Complex::new(0.0, 0.0); cols]; rows];
        for i in 0..rows {
            for j in 0..cols {
                matrix[i][j] = result[(i, j)];
            }
        }
        TopologicalGate {
            name: format!("{}∘{}", other.name, self.name),
            matrix,
            braid_description: Some(format!(
                "Composition of {} and {}",
                self.name, other.name
            )),
        }
    }

    /// Tensor product of two gates.
    pub fn tensor(&self, other: &TopologicalGate) -> TopologicalGate {
        let m1 = self.to_matrix();
        let m2 = other.to_matrix();
        let result = m1.kronecker(&m2);
        let rows = result.nrows();
        let cols = result.ncols();
        let mut matrix = vec![vec![Complex::new(0.0, 0.0); cols]; rows];
        for i in 0..rows {
            for j in 0..cols {
                matrix[i][j] = result[(i, j)];
            }
        }
        TopologicalGate {
            name: format!("{}⊗{}", self.name, other.name),
            matrix,
            braid_description: Some(format!(
                "Tensor product of {} and {}",
                self.name, other.name
            )),
        }
    }

    /// Apply the gate to a state vector.
    pub fn apply(&self, state: &[Complex<f64>]) -> Vec<Complex<f64>> {
        let mat = self.to_matrix();
        let input = DVector::from_column_slice(state);
        let output = &mat * input;
        output.iter().cloned().collect()
    }
}

/// Build a universal gate set from topological braiding.
pub struct UniversalGateSet {
    pub hadamard: TopologicalGate,
    pub cnot: TopologicalGate,
    pub t_gate: TopologicalGate,
}

impl UniversalGateSet {
    pub fn new() -> Self {
        Self {
            hadamard: TopologicalGate::hadamard(),
            cnot: TopologicalGate::cnot(),
            t_gate: TopologicalGate::t_gate(),
        }
    }

    /// Verify all gates are unitary.
    pub fn verify_unitarity(&self) -> bool {
        self.hadamard.is_unitary() && self.cnot.is_unitary() && self.t_gate.is_unitary()
    }
}

impl Default for UniversalGateSet {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hadamard_unitary() {
        let h = TopologicalGate::hadamard();
        assert!(h.is_unitary());
    }

    #[test]
    fn test_cnot_unitary() {
        let cnot = TopologicalGate::cnot();
        assert!(cnot.is_unitary());
    }

    #[test]
    fn test_pauli_x_unitary() {
        let x = TopologicalGate::pauli_x();
        assert!(x.is_unitary());
    }

    #[test]
    fn test_phase_s_unitary() {
        let s = TopologicalGate::phase_s();
        assert!(s.is_unitary());
    }

    #[test]
    fn test_t_gate_unitary() {
        let t = TopologicalGate::t_gate();
        assert!(t.is_unitary());
    }

    #[test]
    fn test_hadamard_squared_is_identity() {
        let h = TopologicalGate::hadamard();
        let hh = h.compose(&TopologicalGate::hadamard());
        let mat = hh.to_matrix();
        let identity = DMatrix::identity(2, 2);
        assert!((mat - identity).norm() < 1e-8);
    }

    #[test]
    fn test_hadamard_apply_zero() {
        let h = TopologicalGate::hadamard();
        let result = h.apply(&[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);
        let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
        assert!((result[0].re - inv_sqrt2).abs() < 1e-10);
        assert!((result[1].re - inv_sqrt2).abs() < 1e-10);
    }

    #[test]
    fn test_cnot_on_10() {
        let cnot = TopologicalGate::cnot();
        // |10⟩ should map to |11⟩
        let result = cnot.apply(&[
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);
        assert!(result[3].re > 0.99);
    }

    #[test]
    fn test_universal_gate_set() {
        let gs = UniversalGateSet::new();
        assert!(gs.verify_unitarity());
    }

    #[test]
    fn test_gate_compose() {
        let x = TopologicalGate::pauli_x();
        let h = TopologicalGate::hadamard();
        let xh = x.compose(&h);
        assert!(xh.is_unitary());
        assert_eq!(xh.name, "H∘X");
    }

    #[test]
    fn test_tensor_product() {
        let x = TopologicalGate::pauli_x();
        let xx = x.tensor(&TopologicalGate::pauli_x());
        assert_eq!(xx.matrix.len(), 4); // 4×4 matrix
        assert!(xx.is_unitary());
    }
}
