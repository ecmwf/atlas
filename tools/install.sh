#! /usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -e -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPT_DIR:$PATH

### COMMAND LINE OPTIONS

# Some defaults for the arguments
PREFIX=$(pwd)/install
CMAKE_OPTIONS=""
with_fckit=false
with_gridtools=false
with_ectrans=false
with_fftw=false
with_qhull=false
with_lz4=false
with_deps=false
with_atlas_orca=false
with_atlas_fesom=false
WORK_DIR=$(pwd)
BUILD_TYPE=RelWithDebInfo

function print_help {
    echo "Quick installer for Atlas and its dependencies"
    echo "----------------------------------------------"
    echo ""
    echo "Usage:"
    echo ""
    echo "  install.sh [--with-deps] [--prefix <prefix>] [--build-type <build-type>] [--cmake <cmake>] [--parallel <nthread>] \\"
    echo "             [--with-atlas-orca] [--with-atlas-fesom] \\"
    echo "             [--enable-fortran] [--enable-gridtools] [--enable-ectrans] [--enable-fftw] [--enable-qhull] [--work-dir <work-dir>] [--help]"
    echo ""
    echo "  "
    echo ""
    echo "Options:"
    echo ""
    echo "  --with-deps                  Install dependencies together with Atlas: fftw, qhull, ecbuild, eckit, fckit, fiat, ectrans"
    echo "  --prefix <prefix>            Install prefix for atlas (and its dependencies if requested with '--with-deps')"
    echo "  --build-type <build-type>    Build type for atlas (and its dependencies if requested with '--with-deps')"
    echo "                               Possible values are ( Release | RelWithDebInfo | Debug )"
    echo "  --cmake <cmake>              Extra CMake Options to configure atlas and its dependencies"
    echo "  --enable-ectrans             Enable optional trans dependency"
    echo "  --enable-gridtools           Enable optional gridtools dependency (only has effect when '--with-deps' is also used)"
    echo "                               ! Requires Boost ! Specify Boost_ROOT environment variable to a fairly recent Boost installation"
    echo "  --enable-fftw                Enable optional fftw dependency required for spectral transforms"
    echo "  --enable-qhull               Enable optional qhull required for meshing of unstructured grids"
    echo "  --enable-fortran             Enable optional fortran interfaces"
    echo "  --with-atlas-orca            Enable optional atlas-orca plugin for ORCA grids"
    echo "  --with-atlas-fesom           Enable optional atlas-fesom plugin for FESOM grids"
    echo "  --work-dir <workdir>         Working directory where sources and builds live"
    echo "  --help                       Print this help"
    echo ""
    echo "Notes:"
    echo ""
    echo "  For each dependency, master branches are used."
    echo "  If different versions are required, please adapt accordingly"
    echo "  Certain dependency features which atlas does not require may be disabled"
    echo ""
    echo "  You may need to set environment variables like 'CC', 'CXX', 'FC' to the desired compilers"
    echo ""
    echo "  "
}

# Parse command line arguments
while [ $# != 0 ]; do
    case "$1" in
    "--cmake")
        CMAKE_OPTIONS="$2"; shift
        ;;
    "--with-deps")
        with_deps=true;
        ;;
    "--enable-fortran")
        with_fckit=true;
        ;;
    "--enable-fftw")
        with_fftw=true;
        ;;
    "--enable-trans")
        with_trans=true;
        ;;
    "--enable-gridtools")
        with_gridtools=true;
        ;;
    "--enable-qhull")
        with_qhull=true;
        ;;
    "--enable-lz4")
        with_lz4=true;
        ;;
    "--with-atlas-orca")
        with_atlas_orca=true;
        ;;
    "--with-atlas-fesom")
        with_atlas_fesom=true;
        ;;
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--build-type")
        BUILD_TYPE="$2"; shift
        ;;
    "--work-dir")
        WORK_DIR="$2"; shift
        ;;
    "--parallel")
        export CMAKE_BUILD_PARALLEL_LEVEL="$2"; shift
        ;;
    "--help")
        print_help;
        exit 0
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

### START OF SCRIPT

echo "Sources and builds will be located in ${WORK_DIR}"

SOURCES_DIR=${WORK_DIR}/sources
BUILDS_DIR=${WORK_DIR}/builds
export PATH=${PREFIX}/bin:${PATH}

export CMAKE_PREFIX_PATH=${PREFIX}:${CMAKE_PREFIX_PATH}

mkdir -p ${SOURCES_DIR}
mkdir -p ${BUILDS_DIR}

if ${with_deps}; then

  ### Install FFTW (optional, off by default)
  if ${with_fftw}; then
    install-fftw.sh --prefix ${PREFIX}
  fi

  ### Install qhull (optional, off by default)
  if ${with_qhull}; then
    install-qhull.sh --prefix ${PREFIX}
  fi

  ### Install lz4 (optional, off by default)
  if ${with_lz4}; then
    install-lz4.sh --prefix ${PREFIX}
  fi

  ### Install ecbuild
  echo "Installing ecbuild"
  [[ -d ${SOURCES_DIR}/ecbuild ]] || git clone -b master https://github.com/ecmwf/ecbuild ${SOURCES_DIR}/ecbuild
  cmake -S ${SOURCES_DIR}/ecbuild -B ${BUILDS_DIR}/ecbuild -DCMAKE_INSTALL_PREFIX=${PREFIX} -DENABLE_TESTS=OFF 
  cmake --build   ${BUILDS_DIR}/ecbuild
  cmake --install ${BUILDS_DIR}/ecbuild

  ### Install eckit
  echo "Installing eckit"
  [[ -d ${SOURCES_DIR}/eckit ]] || git clone -b master https://github.com/ecmwf/eckit ${SOURCES_DIR}/eckit
  cmake -S ${SOURCES_DIR}/eckit -B ${BUILDS_DIR}/eckit \
        -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DENABLE_TESTS=OFF \
        -DENABLE_ECKIT_GEO=OFF \
        -DENABLE_ECKIT_SQL=OFF \
        -DENABLE_ECKIT_CMD=OFF \
        -DENABLE_EIGEN=OFF \
        -DENABLE_LAPACK=OFF \
        -DENABLE_ARMADILLO=OFF \
        -DENABLE_VIENNACL=OFF \
        -DENABLE_CUDA=OFF \
        -DENABLE_AEC=OFF \
        -DENABLE_XXHASH=OFF \
        -DENABLE_JEMALLOC=OFF \
        -DENABLE_BZIP2=OFF \
        -DCMAKE_DISABLE_FIND_PACKAGE_Doxygen=ON \
        ${CMAKE_OPTIONS}
  cmake --build   ${BUILDS_DIR}/eckit
  cmake --install ${BUILDS_DIR}/eckit

  ### Install fckit
  if ${with_fckit}; then
    echo "Installing fckit"
    [[ -d ${SOURCES_DIR}/fckit ]] || git clone -b master https://github.com/ecmwf/fckit ${SOURCES_DIR}/fckit
    cmake ${SOURCES_DIR}/fckit -B ${BUILDS_DIR}/fckit -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_TESTS=OFF \
          ${CMAKE_OPTIONS}
    cmake --build ${BUILDS_DIR}/fckit
    cmake --install ${BUILDS_DIR}/fckit
  fi

  ### Install fiat + ectrans (optional, off by default)
  if ${with_ectrans}; then
    echo "Installing fiat"
    [[ -d ${SOURCES_DIR}/fiat ]] || git clone -b main https://github.com/ecmwf-ifs/fiat ${SOURCES_DIR}/fiat
    cmake -S ${SOURCES_DIR}/fiat -B ${BUILDS_DIR}/fiat \
          -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_TESTS=OFF \
          ${CMAKE_OPTIONS}
    cmake --build   ${BUILDS_DIR}/fiat
    cmake --install ${BUILDS_DIR}/fiat

    echo "Installing ectrans"
    [[ -d ${SOURCES_DIR}/ectrans ]] || git clone -b main https://github.com/ecmwf-ifs/ectrans ${SOURCES_DIR}/ectrans
    cmake -S ${SOURCES_DIR}/ectrans -B ${BUILDS_DIR}/ectrans -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_TESTS=OFF \
          ${CMAKE_OPTIONS}
    cmake --build ${BUILDS_DIR}/ectrans
    cmake --install ${BUILDS_DIR}/ectrans
  fi


  ### Install gridtools (optional, off by default)
  if ${with_gridtools}; then
    echo "Installing gridtools"
    # Note: known to work version: 80187f11
    [[ -d ${SOURCES_DIR}/gridtools ]] || git clone -b master https://github.com/gridtools/gridtools ${SOURCES_DIR}/gridtools
	( cd ${SOURCES_DIR}/gridtools && git checkout 80187f11 )
    cmake -S ${SOURCES_DIR}/gridtools -B ${BUILDS_DIR}/gridtools \
          -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DBUILD_TESTING=OFF \
          -DGT_ENABLE_OPENMP=OFF \
          ${CMAKE_OPTIONS}
    cmake --build   ${BUILDS_DIR}/gridtools
    cmake --install ${BUILDS_DIR}/gridtools
    ### Fix non-standard GridTools installation detection
    if [[ -f ${PREFIX}/lib/cmake/GridToolsConfig.cmake ]]; then
    export GridTools_DIR=${PREFIX}/lib/cmake # see GridTools issue (https://github.com/GridTools/gridtools/issues/1395)
    fi
  fi
fi


### Install atlas
echo "Installing atlas"
cmake -S ${SCRIPT_DIR}/.. -B ${BUILDS_DIR}/atlas \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DENABLE_SANDBOX=ON \
      -DENABLE_TESTS=OFF \
      ${CMAKE_OPTIONS}
cmake --build   ${BUILDS_DIR}/atlas
cmake --install ${BUILDS_DIR}/atlas

if ${with_atlas_orca}; then
  echo "Installing atlas-orca"
  [[ -d ${SOURCES_DIR}/atlas-orca ]] || git clone -b master https://github.com/ecmwf/atlas-orca ${SOURCES_DIR}/atlas-orca
  cmake -S ${SOURCES_DIR}/atlas-orca -B ${BUILDS_DIR}/atlas-orca -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        ${CMAKE_OPTIONS}
  cmake --build   ${BUILDS_DIR}/atlas-orca
  cmake --install ${BUILDS_DIR}/atlas-orca
fi

if ${with_atlas_fesom}; then
  if ${with_deps}; then
    echo "atlas-fesom plugin requires METIS. Installing METIS as dependency."
    install-metis.sh --prefix ${PREFIX}
  fi
  echo "Installing atlas-fesom"
  [[ -d ${SOURCES_DIR}/atlas-fesom ]] || git clone -b ${ATLAS_FESOM_VERSION:-main} ${ATLAS_FESOM_GIT:-"https://github.com/ecmwf/atlas-fesom"} ${SOURCES_DIR}/atlas-fesom
  cmake -S ${SOURCES_DIR}/atlas-fesom -B ${BUILDS_DIR}/atlas-fesom -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        ${CMAKE_OPTIONS}
  cmake --build   ${BUILDS_DIR}/atlas-fesom
  cmake --install ${BUILDS_DIR}/atlas-fesom
fi

