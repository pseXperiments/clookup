# clookup
Lookup argument with function composition

## Benchmark
1. To run benchmark, before running it install `flamegraph`.
```sh
cargo install flamegraph
```

2. In Linux, install `perf`. \
In Mac OS, it uses built-in dtrace.

3. Run benchmark
```sh
make bench
```

4. The result of benchmark will be stored as `flamegraph.svg` in the root.
