#!/bin/bash
cargo build --profile prod \
    && cp target/release/sswsort . \
    || echo "You may need to install Rust nightly!"
