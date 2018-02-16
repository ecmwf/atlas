Atlas
=====

Project home: https://software.ecmwf.int/wiki/display/ATLAS/Atlas
Contact: Willem Deconinck (willem.deconinck@ecmwf.int)
Publication (open-access): https://doi.org/10.1016/j.cpc.2017.07.006

Atlas is a ECMWF framework for parallel data-structures supporting unstructured
grids and function spaces, with the aim to investigate alternative more scalable
dynamical core options for Earth System models, and to support modern interpolation
and product generation softrware

Atlas is predominantly C++ code, with main features available to Fortran codes
through a F2003 interface. It requires some flavour of Unix (such as Linux).
It is known to run on a number of systems, some of which are directly supported
by ECMWF.

If Atlas has been useful in your research, please cite the above-mentioned publication.

Tested compilers include

- GCC 4.9.1, 5.3.0, 6.3.0, 7.2.0
- Intel 15.0.2, 16.0.3, 17.
- CCE 8.4.5, 8.5.8, 8.6.2
- PGI-Fortran 17.7 combined with GNU-C/C++ 6.3

Third Party Requirements
------------------------

Requirements to compile atlas:

- CMake --- For use and installation see http://www.cmake.org/
- ecbuild --- ECMWF library of CMake macros
- eckit (with mpi support)

Recommended:

- fckit --- For enabling Fortran interfaces
- python (only when Fortran bindings are required)

Optional for certain features:

- gridtools-storage --- For GPU interoperability
- transi --- For enabling IFS spherical harmonics transforms ( not open-source )
- CGAL --- For enabling Delaunay triangulation of unstructured grids

Installation
------------

Atlas employs an out-of-source build/install based on CMake.

Make sure ecbuild, eckit and fckit are installed and the ecbuild
executable script is found ( `which ecbuild` ):

```bash
# For finding eckit
ECKIT_PATH               # Path to eckit prefix

# For finding fckit
FCKIT_PATH               # Path to fckit prefix
```

Other environment variables that could be required for optional features:

```bash
# For finding gridtools-storage
GRIDTOOLS_STORAGE_PATH   # Path to gridtools-storage prefix

# For finding transi
TRANSI_PATH              # Path to transi prefix

# For finding CGAL
BOOST_ROOT               # Path to Boost prefix
CGAL_DIR                 # Path to directory containing CGALConfig.cmake
```

Now proceed with installation as follows

```bash
# Environment --- Edit as needed
ATLAS_SRC=$(pwd)
ATLAS_BUILD=build
ATLAS_INSTALL=$HOME/local

# 1. Create the build directory:
mkdir $ATLAS_BUILD
cd $ATLAS_BUILD

# 2. Run CMake
ecbuild --prefix=$ATLAS_INSTALL -- $ATLAS_SRC

# 3. Compile / Install
make -j10
make install

# 4. Check installation
$ATLAS_INSTALL/bin/atlas --info
```

Extra flags maybe added to step 2 to fine-tune configuration.

- `--build=DEBUG|RELEASE|BIT` --- Optimisation level
  * DEBUG:   No optimisation (-O0 -g)
  * BIT:     Maximum optimisation while remaning bit-reproducible (-O2 -g)
  * RELEASE: Maximum optimisation (-O3)
- `-DENABLE_OMP=OFF` --- Disable OpenMP
- `-DENABLE_FORTRAN=OFF` --- Disable Compilation of Fortran bindings

> **Note:**
> By default compilation is done using shared libraries. Some systems have linking
> problems with static libraries that have not been compiled with `-fPIC`.
> In this case, also compile atlas using static linking, by adding to step 2:
    `--static`

Runtime Configuration
---------------------

Atlas behaviour can be configured through some environment variables with defaults marked in square brackets

- `ATLAS_INFO=<0|[1]>`  --- Control printing of Atlas standard information
- `ATLAS_DEBUG=<[0]|1>` --- Control printing of Atlas debug information
- `ATLAS_TRACE=<[0]|1>` --- Control printing of Atlas traces (includes timings)
