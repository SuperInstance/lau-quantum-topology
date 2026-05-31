# lau-quantum-topology

Anyons are particles in 2D whose exchange isn't just a sign flip — it's a phase, and that phase depends on the *path* taken. Braid them, measure the result, and you've performed a computation that's topologically protected. Even a hermit crab scrambling the path can't change the outcome.

This crate implements the mathematics: braid groups, anyon fusion rules, topological quantum gates, Jones polynomials, TQFT, and modular tensor categories.

## The math in 60 seconds

The **braid group** Bₙ has generators σ₁...σₙ₋₁ (swapping adjacent strands) with relations σᵢσⱼ = σⱼσᵢ for |i-j|>1 and σᵢσᵢ₊₁σᵢ = σᵢ₊₁σᵢσᵢ₊₁. These aren't just abstract — they model anyon braiding in 2D topological phases.

Key structures:

- **Fibonacci anyons:** τ × τ = 1 + τ, the simplest universal anyon for TQC
- **Ising anyons:** σ × σ = 1 + ψ, relevant to Majorana zero modes
- **Jones polynomial:** V_L(t) computed from braid representations — a knot invariant from physics
- **TQFT:** a functor from the cobordism category to Vect, assigning vector spaces to surfaces
- **Modular tensor categories:** S-matrix, T-matrix, fusion rules, Verlinde formula

References: Wang, *Topological Quantum Computation* (2010); Kitaev, *Fault-tolerant quantum computation by anyons* (2003)

## Quick start

```rust
use lau_quantum_topology::{
    AnyonType, BraidGroup, TopologicalGate, ModularTensorCategory
};

// Fibonacci anyon fusion: τ ⊗ τ = 1 ⊕ τ
let tau = AnyonType::Fibonacci;
let fusion = tau.fuse(&tau); // {1: 1.0, τ: 1.0}

// Create a braid group B₃
let b3 = BraidGroup::new(3);

// Build a braid: σ₁σ₂⁻¹σ₁
let braid = b3.word(&[1, -2, 1]);

// Burau representation (gives a matrix)
let repr = braid.burau_representation();

// Compute Jones polynomial of the braid closure
let jones = braid.jones_polynomial(-1.0);

// Modular tensor category for Fibonacci anyons
let mtc = ModularTensorCategory::fibonacci();
assert!(mtc.verify_verlinde());    // Verlinde formula
assert!(mtc.verify_modular_ST());  // (ST)³ = S²
```

## Key types

| Type | What it is |
|------|-----------|
| `AnyonType` | Fibonacci or Ising anyon with fusion rules and R/F matrices |
| `BraidGroup` | Bₙ with generators, relations, and representations |
| `TopologicalGate` | Quantum gates from braiding (Hadamard, CNOT, etc.) |
| `JonesPolynomial` | V_L(t) from Burau or Kauffman bracket |
| `TQFT` | 2D and 3D topological quantum field theory functor |
| `ModularTensorCategory` | S/T matrices, fusion rules, quantum dimensions |

## Contributing

[Open an issue](https://github.com/SuperInstance/lau-quantum-topology/issues) or PR. We'd love:

- More anyon models (SU(2)ₖ for k>2)
- Efficient braid simplification algorithms
- Connections to condensed matter simulation
