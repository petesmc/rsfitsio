#!/bin/bash

cd src/bin/testprog

cargo run --bin testprog > rust_testprog.lis

file1="rust_testprog.lis"
file2="testprog.out"

if cmp -s "$file1" "$file2"; then
    printf 'The file "%s" is the same as "%s"\n' "$file1" "$file2"
else
    printf 'The file "%s" is different from "%s"\n' "$file1" "$file2"
    exit 1
fi

file1="testprog.fit"
file2="testprog.std"

if cmp -s "$file1" "$file2"; then
    printf 'The file "%s" is the same as "%s"\n' "$file1" "$file2"
else
    printf 'The file "%s" is different from "%s"\n' "$file1" "$file2"
    exit 1
fi
