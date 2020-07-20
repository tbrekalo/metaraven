# Metaraven

Metaraven is a de novo genome assembler for long uncorrected reads 
on metagenomic populations. 

## Usage
To build metaraven run the following commands:
```bash
git clone --recursive https://github.com/tbrekalo/metaraven.git metaraven
cd metaraven && bash ./build.sh

```
which will display the following usage:
```bash
usage: raven [options ...] <sequences>

# default output is stdout
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    --weaken
      use larger (k, w) when assembling highly accurate sequences
    -p, --polishing-rounds <int>
      default: 2
      number of times racon is invoked
    -m, --match <int>
      default: 3
      score for matching bases
    -n, --mismatch <int>
      default: -5
      score for mismatching bases
    -g, --gap <int>
      default: -4
      gap penalty (must be negative)
    --graphical-fragment-assembly <string>
      prints the assemblg graph in GFA format
    --resume
      resume previous run from last checkpoint
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
       prints the usage

  only available when built with CUDA:
    -c, --cuda-poa-batches <int>
       default: 0
       number of batches for CUDA accelerated polishing
    -b, --cuda-banded-alignment
       use banding approximation for polishing on GPU
       (only applicable when -c is used)
    -a, --cuda-alignment-batches <int>
       default: 0
       number of batches for CUDA accelerated alignment
```

#### Dependencies
- gcc 4.8+ or clang 3.5+
- cmake 3.9+
- zlib

### CUDA Support
To build submodule racon with CUDA support, add `-Dracon_enable_cuda=ON` while running `cmake`. For more information see [this](https://github.com/lbcb-sci/racon).

#### Dependnecies
- gcc 5.0+
- cmake 3.10+
- CUDA 9.0+

## Acknowledgment

This work has been supported in part by the Croatian Science Foundation under the projects Algorithms for genome sequence analysis (UIP-11-2013-7353) and Single genome and metagenome assembly (IP-2018-01-5886), and in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS).
