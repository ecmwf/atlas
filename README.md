Atlas
=====

Project home: https://software.ecmwf.int/wiki/display/ATLAS/Atlas
Contact: Willem Deconinck (willem.deconinck@ecmwf.int)

Atlas is a ECMWF framework for parallel data-structures supporting unstructured
grids and function spaces, with the aim to investigate alternative more scalable
dynamical core options for IFS, and to create modern interpolation and product
generation softrware

Atlas is predominantly C++ code, with main features available to Fortran codes
through a F2003 interface. It requires some flavour of Unix (such as Linux).
It is known to run on a number of systems, some of which are directly supported
by ECMWF.

Tested compilers include

- GCC 4.8.1
- Intel 13.0.1, 14.0.1
- CCE 8.2.7, 8.3.1


Third Party Requirements
------------------------

Requirements to compile atlas:

- CMake --- For use and installation see http://www.cmake.org/
- ecbuild
- eckit (with mpi support)
- python (only when Fortran bindings are required)

Recommended:

- OpenMP --- Shared memory parallelisation
- fkit --- Unit testing for Fortran


Installation
------------

Atlas employs an out-of-source build/install based on CMake.

```bash
# Environment --- Edit as needed
ATLAS_SRC=$(pwd)
ATLAS_BUILD=build
ATLAS_INSTALL=$HOME/local
export CC=mpicc
export CXX=mpicxx

# 1. Create the build directory:
mkdir $ATLAS_BUILD
cd $ATLAS_BUILD

# 2. Run CMake
cmake $ATLAS_SRC -DCMAKE_INSTALL_PREFIX=$ATLAS_INSTALL

# 3. Compile / Install
make -j10
make install
```

Extra flags maybe added to step 2 to fine-tune configuration.

- `-DCMAKE_BUILD_TYPE=DEBUG|RELEASE|BIT` --- Optimisation level
  * DEBUG:   No optimisation
  * BIT:     Maximum optimisation while remaning bit-reproducible
  * RELEASE: Maximum optimisation
- `-DENABLE_OMP=OFF` --- Disable OpenMP
- `-DENABLE_MPI=OFF` --- Disable MPI
- `-DENABLE_FORTRAN=OFF` --- Disable Compilation of Fortran bindings

> **Note:**
> By default compilation is done using shared libraries. Some systems have linking
> problems with static libraries that have not been compiled with `-fPIC`.
> In this case, also compile atlas using static linking, by adding to step 2:
    `-DBUILD_SHARED_LIBS=OFF`

Runtime Configuration
---------------------

Atlas behaviour can be configured through variables defined at the command-line, in the
environment, or in configuration files.
In following table, the column "variable" can be edited in configuration files.

| variable                    | command line     | environment              | default            |
|-----------------------------|------------------|--------------------------|--------------------|
| `debug`                     | `--debug`        | `$DEBUG`                 | `0`                |
|                             | `--name`         |                          | `<app>`            |
| `atlas.logfile`             | `--logfile`      | `$ATLAS_LOGFILE`         |                    |
| `atlas.logfile_task`        | `--logfile_task` | `$ATLAS_LOGFILE_TASK`    | `-1`               |
| `atlas.console_task`        | `--console_task` | `$ATLAS_CONSOLE_TASK`    | `0`                |
