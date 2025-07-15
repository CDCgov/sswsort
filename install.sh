#!/bin/bash
cargo build --release \
    && cp target/release/sswsort . \
    || echo "You may need to install Rust nightly!"
