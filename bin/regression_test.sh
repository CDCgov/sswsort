#!/bin/bash

VERSION="v1.3"

for tag in "$VERSION" master; do
    echo "Running '$tag'" \
        && git checkout $tag > /dev/null \
        && time ./SSWSORT flu-ABCD90P presets/flu-ABCD90P.fasta test-flu-$tag.txt \
        && time ./SSWSORT cov-beta presets/cov-beta.fasta test-cov-$tag.txt
done

for mod in flu cov; do
    echo "Checking '$mod'" \
        && diff <(cut -f2- test-${mod}-master.txt | sort) <(cut -f2- test-${mod}-${VERSION}.txt | sort) \
        && echo "  Passed '$mod'" \
        || echo "  Failed '$mod'"
done
