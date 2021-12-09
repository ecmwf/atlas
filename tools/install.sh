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
with_gridtools=false
with_trans=false
with_fftw=false
with_cgal=false
with_deps=false
WORK_DIR=$(pwd)
BUILD_TYPE=RelWithDebInfo

function print_help {
    echo "Quick installer for Atlas and its dependencies"
    echo "----------------------------------------------"
    echo ""
    echo "Usage:"
    echo ""
    echo "  install.sh [--with-deps] [--prefix <prefix>] [--build-type <build-type>] [--cmake <cmake>] \\"
    echo "             [--enable-gridtools] [--enable-trans] [--enable-fftw] [--enable-cgal] [--work-dir <work-dir>] [--help]"
    echo ""
    echo "  "
    echo ""
    echo "Options:"
    echo ""
    echo "  --with-deps                  Install dependencies together with Atlas: fftw, cgal, ecbuild, eckit, fckit, fiat, ectrans"
    echo "  --prefix <prefix>            Install prefix for atlas (and its dependencies if requested with '--with-deps')"
    echo "  --build-type <build-type>    Build type for atlas (and its dependencies if requested with '--with-deps')"
    echo "                               Possible values are ( Release | RelWithDebInfo | Debug )"
    echo "  --cmake <cmake>              Extra CMake Options to configure atlas and its dependencies"
    echo "  --enable-trans               Enable optional trans dependency"
    echo "  --enable-gridtools           Enable optional gridtools dependency (only has effect when '--with-deps' is also used)"
    echo "                               ! Requires Boost ! Specify Boost_ROOT environment variable to a fairly recent Boost installation"
    echo "  --enable-fftw                Enable optional fftw dependency required for spectral transforms"
    echo "  --enable-cgal                Enable optional cgal required for meshing of unstructured grids"
    echo "                               ! Requires Boost ! Specify Boost_ROOT environment variable to a fairly recent Boost installation"
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
    "--enable-fftw")
        with_fftw=true;
        ;;
    "--enable-trans")
        with_trans=true;
        ;;
    "--enable-gridtools")
        with_gridtools=true;
        ;;
    "--enable-cgal")
        with_cgal=true;
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

export CMAKE_PREFIX_PATH=${PREFIX}

mkdir -p ${SOURCES_DIR}
mkdir -p ${BUILDS_DIR}

if ${with_deps}; then

  ### Install FFTW (optional, off by default)
  if ${with_fftw}; then
    install-fftw.sh --prefix ${PREFIX}
  fi

  ### Install CGAL (optional, off by default)
  if ${with_cgal}; then
    echo "Installing CGAL"
    [[ -d ${SOURCES_DIR}/cgal ]] || git clone https://github.com/CGAL/cgal.git ${SOURCES_DIR}/cgal
    cd ${SOURCES_DIR}/cgal
    git checkout releases/CGAL-5.0
    mkdir -p ${BUILDS_DIR}/cgal && cd ${BUILDS_DIR}/cgal
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_TESTING=OFF \
          ${SOURCES_DIR}/cgal
    make -j8 install
  fi

  ### Install ecbuild
  echo "Installing ecbuild"
  [[ -d ${SOURCES_DIR}/ecbuild ]] || git clone -b master https://github.com/ecmwf/ecbuild ${SOURCES_DIR}/ecbuild
  mkdir -p ${BUILDS_DIR}/ecbuild && cd ${BUILDS_DIR}/ecbuild
  cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DENABLE_TESTS=OFF ${SOURCES_DIR}/ecbuild
  make -j8 install

  ### Install eckit
  echo "Installing eckit"
  [[ -d ${SOURCES_DIR}/eckit ]] || git clone -b master https://github.com/ecmwf/eckit ${SOURCES_DIR}/eckit
  mkdir -p ${BUILDS_DIR}/eckit && cd ${BUILDS_DIR}/eckit
  cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DENABLE_TESTS=OFF \
        -DENABLE_ECKIT_SQL=OFF \
        -DENABLE_ECKIT_CMD=OFF \
        -DENABLE_EIGEN=OFF \
        -DENABLE_LAPACK=OFF \
        -DENABLE_ARMADILLO=OFF \
        -DENABLE_VIENNACL=OFF \
        -DENABLE_CUDA=OFF \
        -DENABLE_AEC=OFF \
        -DENABLE_XXHASH=OFF \
        -DENABLE_LZ4=OFF \
        -DENABLE_JEMALLOC=OFF \
        -DENABLE_BZIP2=OFF \
        -DCMAKE_DISABLE_FIND_PACKAGE_Doxygen=ON \
        ${CMAKE_OPTIONS} \
        ${SOURCES_DIR}/eckit
  make -j8 install

  ### Install fckit
  echo "Installing fckit"
  [[ -d ${SOURCES_DIR}/fckit ]] || git clone -b master https://github.com/ecmwf/fckit ${SOURCES_DIR}/fckit
  mkdir -p ${BUILDS_DIR}/fckit && cd ${BUILDS_DIR}/fckit
  cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DENABLE_TESTS=OFF \
        ${CMAKE_OPTIONS} \
        ${SOURCES_DIR}/fckit
  make -j8 install

  ### Install fiat + ectrans (optional, off by default)
  if ${with_trans}; then
    echo "Installing fiat"
    [[ -d ${SOURCES_DIR}/fiat ]] || git clone -b master https://github.com/ecmwf-ifs/fiat ${SOURCES_DIR}/fiat
    mkdir -p ${BUILDS_DIR}/fiat && cd ${BUILDS_DIR}/fiat
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_TESTS=OFF \
          ${CMAKE_OPTIONS} \
          ${SOURCES_DIR}/fiat
    make -j8 install

    echo "Installing ectrans"
    [[ -d ${SOURCES_DIR}/ectrans ]] || git clone -b master https://github.com/ecmwf-ifs/ectrans ${SOURCES_DIR}/ectrans
    mkdir -p ${BUILDS_DIR}/ectrans && cd ${BUILDS_DIR}/ectrans
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_TESTS=OFF \
          ${CMAKE_OPTIONS} \
          ${SOURCES_DIR}/ectrans
    make -j8 install
  fi


  ### Install gridtools (optional, off by default)
  if ${with_gridtools}; then
    echo "Installing gridtools"
    # Note: known to work version: 80187f11
    [[ -d ${SOURCES_DIR}/gridtools ]] || git clone -b master https://github.com/gridtools/gridtools ${SOURCES_DIR}/gridtools
	( cd ${SOURCES_DIR}/gridtools && git checkout 80187f11 )
    mkdir -p ${BUILDS_DIR}/gridtools && cd ${BUILDS_DIR}/gridtools
    cmake ${ECBUILD_MODULE_PATH} \
          -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DBUILD_TESTING=OFF \
          -DGT_ENABLE_OPENMP=OFF \
          ${CMAKE_OPTIONS} \
          ${SOURCES_DIR}/gridtools
    make -j8 install
    ### Fix non-standard GridTools installation detection
    if [[ -f ${PREFIX}/lib/cmake/GridToolsConfig.cmake ]]; then
    export GridTools_DIR=${PREFIX}/lib/cmake # see GridTools issue (https://github.com/GridTools/gridtools/issues/1395)
    fi
  fi
fi


### Install atlas
echo "Installing atlas"
mkdir -p ${BUILDS_DIR}/atlas && cd ${BUILDS_DIR}/atlas
cmake ${ECBUILD_MODULE_PATH} \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DENABLE_SANDBOX=ON \
      ${CMAKE_OPTIONS} \
      ${SCRIPT_DIR}/..
make -j8
make install
