# ARTED: Ab-initio Real-Time Electron Dynamics Simulator

- [Overview](#overview)
- [Build](#build)
- [Execution](#execution)
- [Test Environment](#test-environment)
- [License](#license)
- [Acknowledgement](#acknowledgement)

## Overview

ARTED (Ab-initio Real-Time Electron Dynamics simulator) is an open-source
computer program for first-principles calculations of electron dynamics and
light-matter interactions [1,2]. It is based on time-dependent density functional theory
solving time-dependent Kohn-Sham equation in real time using pseudopotentials
and real-space grid representation.

ARTED has been developed in such a way that it runs
optimally in the following supercomputer platforms:

- K-computer [3]
- Fujitsu FX100 supercomputer system [4]
- Linux PC Cluster with x86-64 CPU
- Linux PC Cluster with Intel Knights Landing [4]
- Linux PC Cluster with Intel Knights Corner [5]
- Linux PC Cluster with NVIDIA GPU (OpenACC, Kepler and newer GPUs)

ARTED has been developed by ARTED developers with support from
Center for Computational Sciences, University of Tsukuba.

### Reference
1. G. F. Bertsch, J.-I. Iwata, Angel Rubio, and K. Yabana: "Real-space, real-time method for the dielectric function", [Phys. Rev. B 62, 7998 (2000)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.62.7998).
2. K. Yabana, T. Sugiyama, Y. Shinohara, T. Otobe, and G. F. Bertsch: "Time-dependent density functional theory for strong electromagnetic fields in crystalline solids", [Phys. Rev. B  85, 045134 (2012)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.045134).
3. Shunsuke A. Sato, and Kazuhiro Yabana: "Maxwell + TDDFT multi-scale simulation for laser-matter interaction", [J. Adv. Simulat. Sci. Eng. 1, 98 (2014)](https://www.jstage.jst.go.jp/article/jasse/1/1/1_98/_article).
4. Yuta Hirokawa: "Electron Dynamics Simulation with Time-Dependent Density Functional Theory on Large Scale Many-Core Systems", [SC16 ACM SRC Poster (2016)](http://sc16.supercomputing.org/presentation/?id=spost147&sess=sess408).
5. Yuta Hirokawa, Taisuke Boku, Shunsuke A. Sato, and Kazuhiro Yabana: "Electron Dynamics Simulation with Time-Dependent Density Functional Theory on Large Scale Symmetric Mode Xeon Phi Cluster", [IEEE IPDPS Workshop PDSEC'16 (2016)](http://ieeexplore.ieee.org/abstract/document/7530004/).


## Build

We use [CMake](https://cmake.org/) cross-platform build tools.

CMake detects the following configurations automatically,

- MPI Fortran/C compiler
- OpenMP compile flag
- LAPACK(/BLAS) libraries

CMake software version **must** be 2.8 or later.
We recommend *3.0* or later versions.

### for your computer

    $ mkdir ./build_temp
    $ cd build_temp
    $ ../configure.py && make
    or
    $ mkdir build_temp
    $ cd build_temp
    $ cmake .. && make

You can specify compilers and/or compiler options.

    $ ../configure.py FC=mpiifort CC=mpiicc FFLAGS="-xAVX" CFLAGS="-restrict -xAVX"

### for Supercomputer systems and clusters

We provide the build configuration of the systems with CMake cross-compile mode.

    $ mkdir build_temp
    $ cd build_temp
    $ ../configure.py --arch=<COMPILER>-<SYSTEM> && make
    or
    $ mkdir build_temp
    $ cd build_temp
    $ cmake -D CMAKE_TOOLCHAIN_FILE=<COMPILER>-<SYSTEM> .. && make

For example, the following command builds the application for K-computer.

    $ ../configure.py --arch=fujitsu-k

CMake searches the cross-compile configuration files below `platform` directory.
If you want execution at the system that configuration file is not provided, you can create it yourself.

### Select the simulation mode

ARTED provides two simulation modes, Single-cell (sc) and Multi-scale (ms).
The default target is sc mode. You can select a target with `-t TARGET, --target=TARGET` option in `configure.py` script.
For help type `./configure.py -h`.

    $ ../configure.py -t sc        # build for SC (default)
    $ ../configure.py --target=ms  # build for MS


## Execution

Please read to jobscript directory files.

    $ mpirun -np $NUM_MPI_PROCS ./bin/ARTED_sc.cpu < ./data/input_sc.dat


## Test Environment

### Intel Knights Landing

Oakforest-PACS at JCAHPC, The University of Tokyo and University of Tsukuba

1. Intel Compiler version 17.0.1
2. Intel MPI 5.1.3
3. Intel MKL 11.3.2

### Intel x86-64 CPU with Intel Knights Corner

COMA at CCS, University of Tsukuba

1. Intel Compiler version 16.0.2
1. Intel MPI 5.1.3
1. Intel MKL 11.3.2

### x86-64 CPUs with NVIDIA Kepler GPU

HA-PACS/TCA at CCS, University of Tsukuba

1. PGI Compiler 16.4
2. OpenMPI 1.10.3 or MVAPICH2 GDR 2.1
3. CUDA 7.5.18

### x86-64 CPUs (General version)

1. GCC version 4.4.7
2. OpenMPI 1.10.3
3. LAPACK 3.6.0

### Supercomputer system

K-computer at RIKEN AICS

1. Fujitsu Compiler version K-1.2.0-20-1

FX100 system at Nagoya University

1. Fujitsu Compiler Driver Version 2.0.0


## License

ARTED is available under Apache License version 2.0.

    Copyright 2016 ARTED developers
    
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
    
       http://www.apache.org/licenses/LICENSE-2.0
    
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

## Acknowledgement

### NVIDIA GPU support with OpenACC

Thanks to Mr. Akira Naruse (NVIDIA Corporation)
