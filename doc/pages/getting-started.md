Getting started    {#getting-started}
===============

@tableofcontents

Requirements  {#getting-started-requirements}
------------

Tested compilers include:

- GCC 4.9.1, 5.3.0, 6.3.0, 7.2.0
- Intel 15.0.2, 16.0.3, 17, 18
- CCE 8.4.5, 8.5.8, 8.6.2
- PGI-Fortran 17.7 combined with GNU-C/C++ 6.3
- PGI 17.7

Known compilers to fail include:

- PGI-Fortran 17.10, 18.1

Required dependencies:

- CMake --- For use and installation see http://www.cmake.org/
- ecbuild --- ECMWF library of CMake macros
- eckit (with MPI support) --- C++ support library

Recommended dependencies:

- fckit --- For enabling Fortran interfaces
- python (only when Fortran bindings are required)

Optional dependencies:

- gridtools-storage --- For GPU interoperability
- transi --- For enabling IFS spherical harmonics transforms ( not open-source )
- CGAL --- For enabling Delaunay triangulation of unstructured grids
- Eigen3 -- For certain linear algebra operations
- FFTW -- For enabling inverse spherical harmonics transforms (TransLocal)

Installation  {#getting-started-installation}
------------

Atlas employs an out-of-source build/install based on CMake.

Make sure ecbuild, eckit and fckit are installed and the ecbuild
executable script is found ( `which ecbuild` ). Following environment variables
help the build system to detect the right dependencies:

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
Eigen3_DIR               # Path to directory containing Eigen3Config.cmake
FFTW_PATH                # Path to FFTW prefix
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
- `--build=DEBUG|RELEASE|BIT` --- Optimisation level <dl class="m-diary">
    <dt> DEBUG   </dt>  <dd> No optimisation (`-O0 -g`) </dd>
    <dt> BIT     </dt>  <dd> Maximum optimisation while remaning bit-reproducible (`-O2 -g`) </dd>
    <dt> RELEASE </dt>  <dd> Maximum optimisation (`-O3`) </dd></dl>
- `-DENABLE_OMP=OFF` --- Disable OpenMP
- `-DENABLE_FORTRAN=OFF` --- Disable Compilation of Fortran bindings

@note
    By default compilation is done using shared libraries. Some systems have linking
    problems with static libraries that have not been compiled with `-fPIC`.
    In this case, also compile atlas using static linking, by adding to step 2:
    `--static`

Runtime Configuration   {#getting-started-runtime-configuration}
---------------------

Atlas behaviour can be configured through some environment variables with defaults marked in square brackets

- `ATLAS_INFO=<0|[1]>`  --- Control printing of Atlas standard information
- `ATLAS_DEBUG=<[0]|1>` --- Control printing of Atlas debug information
- `ATLAS_TRACE=<[0]|1>` --- Control printing of Atlas traces (includes timings)


Additional information   {#getting-started-more}
----------------------

-   @subpage building --- @copybrief building
-   @subpage cmake --- @copybrief cmake
-   @subpage custom-buildsystems --- @copybrief custom-buildsystems
