Minimizer benchmark
====

A Benchmark for minimizer evaluation

## Build
Please make sure the following systems have been installed before building from the source.
- GNU/Linux GCC ≥ 4.9.0
- CMAKE ≥ 3.0.0
- zlib ≥ 1.2

```bash
#To build from source, please type in the commandline
mkdir -p build/release && cd $_
CMake [path to source]
make minibench -j4 #use 4 threads to compile
```


