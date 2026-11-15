# LIVIA Rust Core & WebAssembly Engine (`rust/`)

Unified Rust workspace powering LIVIA and LIVIA Atlas across native CLI and browser WebAssembly:

- **`crates/livia-core`**: Shared interface-scoring engine (`iLIS`, `iLIA`, `iLISA`, `LIS`, `cLIS`, `ipSAE`, `actifpTM*`, `pDockQ`/`LIpDockQ`, `pDockQ2`/`LIpDockQ2`, `LIR`/`cLIR`), bitset-accelerated cLIP cosine-distance hierarchical clustering and silhouette selection, Leiden community detection (modularity & CPM), greedy Steiner subnetwork extraction, and deterministic graph layouts (Fruchterman-Reingold, Kamada-Kawai, DrL).
- **`crates/lis-rs`**: Single-binary, multi-core command-line scorer (`lis-rs`) compatible with `python/lis.py` across AlphaFold-Multimer/ColabFold/AlphaPulldown, AlphaFold 3, Boltz-1/2, Chai-1, OpenFold3, Protenix-v2, ESMFold2, and generic coordinate + PAE inputs.
- **`crates/livia-wasm`**: `wasm32-unknown-unknown` target bundled by `scripts/build_wasm.py` into self-contained `../js/livia-wasm.js` (with automatic pure-JavaScript fallback when WebAssembly is blocked by strict CSP).
- **`crates/livia-atlas-cli`**: Streaming `.tar.xz` merger, deduplicator, and `LVTP2` (`tested.bin`) binary index builder for LIVIA Atlas.

## Build & Test

```bash
# Run unit tests across all crates
cargo test --release

# Build native CLI binaries (lis-rs and livia-atlas-cli)
cargo build --release

# Rebuild WebAssembly bundle (writes dist/livia-wasm.js)
cargo build --profile wasm-release --target wasm32-unknown-unknown -p livia-wasm
python3 scripts/build_wasm.py
cp dist/livia-wasm.js ../js/livia-wasm.js
```
