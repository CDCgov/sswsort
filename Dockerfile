FROM redhat/ubi8:latest AS builder

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH

RUN yum update -y && yum install -y zip git which gcc && yum clean all

RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "aarch64" ]; then  RUSTUP_SHA256="9732d6c5e2a098d3521fca8145d826ae0aaa067ef2385ead08e6feac88fa5792"; \
    elif [ "$ARCH" = "x86_64" ]; then RUSTUP_SHA256="4acc9acc76d5079515b46346a485974457b5a79893cfb01112423c89aeb5aa10"; \
    else echo "Unsupported architecture: $ARCH" && exit 1; fi && \
    RUSTUP_URL="https://static.rust-lang.org/rustup/archive/1.29.0/${ARCH}-unknown-linux-gnu/rustup-init" && \
    curl --proto '=https' --tlsv1.2 -sSf -o rustup-init "$RUSTUP_URL" && \
    echo "${RUSTUP_SHA256} *rustup-init" | sha256sum -c - && \
    chmod +x rustup-init && \
    ./rustup-init -y --no-modify-path --profile minimal --default-toolchain nightly && \
    rm rustup-init && \
    chmod -R a+w $RUSTUP_HOME $CARGO_HOME && rustc --version

SHELL ["/bin/bash", "-c"]
WORKDIR /sswsort
ARG sswsort_branch

COPY . .

RUN latest=$(git tag | tail -n1) \
    && git checkout ${sswsort_branch:-$latest} \
    && cargo build --workspace --profile prod \
    && cargo test --workspace

FROM debian:bookworm-slim AS base

ARG APT_MIRROR_NAME=
RUN if [ -n "$APT_MIRROR_NAME" ]; then sed -i.bak -E '/security/! s^https?://.+?/(debian|ubuntu)^http://'"$APT_MIRROR_NAME"'/\1^' /etc/apt/sources.list && grep '^deb' /etc/apt/sources.list; fi
RUN apt-get update --allow-releaseinfo-change --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y procps \
    && apt clean autoclean \
    && apt autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

WORKDIR /app
COPY  --from=builder /sswsort/docs /app/docs
COPY   --from=builder \
    /sswsort/target/prod/sswsort-cli \
    /sswsort/Cargo.toml \
    /sswsort/Cargo.lock \
    /sswsort/LICENSE \
    /sswsort/CHANGELOG.md \
    /sswsort/CITATION.bib \
    /sswsort/CONTRIBUTORS.md \
    /sswsort/README.md \
    /app/

ENV PATH="/app:${PATH}"
WORKDIR /data
