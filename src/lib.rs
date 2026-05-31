#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]
//! # lau-quantum-topology
//!
//! Topological quantum computing primitives for agent reasoning.
//!
//! This crate implements:
//! - Anyon braiding (Fibonacci anyons, Ising anyons)
//! - Braid group representation (B_n generators and relations)
//! - Topological quantum gates from braids (CNOT, Hadamard via braiding)
//! - Jones polynomial computation from braids
//! - TQFT axioms (2D and 3D TQFT functors)
//! - Modular tensor categories (simple objects, fusion rules, S-matrix, T-matrix)
//! - Quantum dimension from fusion rules
//! - Application: agent reasoning as topological quantum computation

pub mod anyon;
pub mod braid;
pub mod gates;
pub mod jones;
pub mod tqft;
pub mod mtc;

pub use anyon::{AnyonType, Anyon, FusionResult};
pub use braid::{Braid, BraidWord};
pub use gates::TopologicalGate;
pub use jones::jones_polynomial;
pub use tqft::{TQFT2D, TQFT3D};
pub use mtc::{ModularTensorCategory, FusionRule, SimpleObject};
