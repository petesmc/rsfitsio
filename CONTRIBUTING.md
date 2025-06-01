# Contributing to rsfitsio

`rsfitsio` welcomes contributions from everyone

## PR Checking

1. Make sure tests pass
2. Run `cargo fmt`
3. Run `cargo fix`
4. Fix new warnings
5. Keep PRs small, and if required split into multiple explainable commits
6. Add tests

## Aims

Initially we are trying to stay close to the `cfitsio` implementation so that its easier to sync fixes. 

This means usage the same function names, variable names, parameter order, etc. In time this will be relaxed.

## Testing

Testing is still very limited and automated testing vs the C-API still needs to be added. In the meantime, it is suggested to compile and run testprog.c and verify that the output matches.

Both diffs should return not differences at all.

```
export LD_LIBRARY_PATH=~/rsfitsio/target/debug

gcc -L ../target/debug -o testprog testprog.c -lrsfitsio

./testprog > out.log

diff out.log testprog.out

diff testprog.fit testprog.std
```

## Sync'ing upstream

Updated `SYNCED_COMMIT.md` file with the commit hash of the upstream cfitsio library when you are syncing.

# Tips

## Testprog

```
cargo run --bin testprog > rust_testprog.lis
```

## Windows

```Powershell
$env:Path += ';C:\code\rsfitsio\target\debug'
```

## Miri

```
RUST_BACKTRACE=1 MIRIFLAGS="-Zmiri-env-forward=RUST_BACKTRACE -Zmiri-disable-isolation -Zmiri-backtrace=full" cargo miri test -- -- tests::test_write_image

MIRIFLAGS="-Zmiri-disable-isolation" cargo miri test

RUST_BACKTRACE=1 MIRIFLAGS="-Zmiri-env-forward=RUST_BACKTRACE -Zmiri-disable-isolation -Zmiri-backtrace=full" cargo miri run --bin testprog
```