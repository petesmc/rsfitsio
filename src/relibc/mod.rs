pub mod header;
// Vendored from the Rust project's std::io (see the copyright headers in
// io/mod.rs and io/buffered.rs). Kept byte-identical to upstream so it can be
// re-synced, so its lints are allowed here rather than fixed in place.
#[allow(clippy::manual_clear, clippy::redundant_guards)]
pub mod io;
pub mod platform;
