# lau-quantum-topology

Topological quantum computing for agent reasoning — anyons, braids, TQFT, and modular tensor categories.

## Overview

This crate implements topological quantum computing primitives that can be used for modeling agent reasoning as topological quantum computation. It provides:

- **Anyon braiding**: Fibonacci anyons and Ising anyons with fusion rules, R-matrices, and F-matrices
- **Braid group representation**: B_n generators σ_i with the braid group relations (Yang-Baxter equation, far commutativity) via the Burau representation
- **Topological quantum gates**: CNOT, Hadamard, Pauli-X, Phase, and T-gates constructed from braiding
- **Jones polynomial**: Computation from braids using the Burau representation
- **TQFT axioms**: 2D and 3D topological quantum field theory functors (Frobenius algebras, Verlinde dimensions)
- **Modular tensor categories**: Simple objects, fusion rules, S-matrix, T-matrix, quantum dimensions, Verlinde formula
- **Agent reasoning**: Braided protocols for modeling topological quantum computation in agent systems

## Usage

```rust
use lau_quantum_topology::*;

// Work with Fibonacci anyons
let result = anyon::fuse_fibonacci(&AnyonType::Fibonacci, &AnyonType::Fibonacci);
// τ × τ = 1 + τ

// Create braids and verify relations
let mut bw = braid::BraidWord::identity(4);
bw.sigma(0).sigma(1).sigma(0);
let b = braid::Braid::new(bw);

// Compute Jones polynomial
let v = jones::jones_trefoil(Complex::new(0.5, 0.0));

// Modular tensor categories
let fib = mtc::ModularTensorCategory::fibonacci();
assert!(fib.verify_s_unitarity());
assert!(fib.verify_verlinde());
```

## Testing

```bash
cargo test
```

70 tests covering braid relations, Jones polynomial values, fusion rules, S-matrix unitarity, Verlinde formula, and modular relations.

## License

MIT
