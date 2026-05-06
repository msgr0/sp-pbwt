# Sampled and Parallel Positional Burrows-Wheeler Transform

## Install

```bash
export LIBOMP=/path/to/lib
export HTSLIB=/path/to/lib
git clone https://github.com/msgr0/sp-pbwt.git
cd sp-pbwt
make
```

## BCF input:
Run `./sp-pbwt-bcf MODE FILE`, with one of the available modes:

**Syscall**:
- linear
- sampled
- blockpar
- stagpar

## BM input:

Run `./sp-pbwt-bm MODE FILE`, with one of the available modes:

**Syscall**:
- linear-syscall
- sampled-syscall
- blockpar-syscall
- stagpar-syscall

**Memory Mapping**:
- linear-mmap
- sampled-mmap
- blockpar-mmap
- stagpar-mmap

