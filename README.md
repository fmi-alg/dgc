# DGC - Distance Closures: Unifying Search- and Lookup-based Shortest Path Speedup Techniques

This Repository holds the source code discussed in our paper "Distance Closures: Unifying Search- and Lookup-based Shortest Path Speedup Techniques".

## Clone

Clone with recursive:

```bash
git clone --recursive https://github.com/fmi-alg/dgc.git
```

## Build

```bash
cmake -B build -S ./cpp -DCMAKE_BUILD_TYPE=Release -DOPTIMIZE_ULTRA=TRUE
cmake --build build
```

### Building with support for no-mmap and direct-io

```bash
cmake -B build -S ./cpp -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE -DSSERIALIZE_INLINE_IN_LTO_ENABLED=TRUE -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Note: gcc-10/g++-10 is needed to compile sserialize.

## Run

Sample data can be obtained from [here](http://data.oscar-web.de/dgc).

### In-memory

```bash
build/dgc boost -i <file>
```

### Export to offset array variant

```bash
build/dgc export-oa -i <file>
```

### Run offset array variant

```bash
build/dgc oa -i data.sserialize-oa
```

### Poor man's async

```bash
taskset 0x1 ./build/dgc -j 8 ...
```

This pins the process to core 1 (hence a single "real" thread) whereas 8 threads are spawned.
If a thread waits for io then the OS will switch to an other one.
