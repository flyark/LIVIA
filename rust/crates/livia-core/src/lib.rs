#![forbid(unsafe_code)]

//! Unified pure-Rust computational core for LIVIA (`lis-rs` CLI & browser Wasm)
//! and LIVIA Atlas (`clipworker.js`, `app.js` network builder, offline indexers).
//!
//! Zero external dependencies so this crate compiles cleanly to both native
//! multi-core targets (`x86_64`/`aarch64` macOS, Linux `gnu`/`musl`, Windows)
//! and tiny `wasm32-unknown-unknown` binaries without C toolchains.

pub mod clip;
pub mod graph;
pub mod metrics;
pub mod structure;
