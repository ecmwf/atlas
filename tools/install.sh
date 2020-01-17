#! /usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -e -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

# Some defaults for the arguments
PREFIX=$(pwd)/install
CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
with_gridtools=false
with_deps=false
if [ -z "${TMPDIR+x}" ]; then
  TMPDIR=${HOME}/tmp
fi
cleanup=""

function print_help {
    echo "install.sh [--cmake <options>] [--with-deps] [--enable-gridtools] [--prefix <prefix>] [--cleanup] [--help]"
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
    "--enable-gridtools")
        with_gridtools=true;
        ;;
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--tmpdir")
        TMPDIR="$2"; shift
        ;;
    "--cleanup")
        cleanup="--cleanup";
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

echo "Downloads and builds will be located in ${TMPDIR}"

if ${with_deps}; then

  ### Install FFTW
  install-fftw.sh --prefix ${PREFIX}

  ### Install ecbuild
  install-dep.sh --repo ecbuild --prefix ${PREFIX} ${cleanup} --cmake "${CMAKE_OPTIONS}"
  export ECBUILD_MODULE_PATH=${PREFIX}/share/ecbuild/cmake

  ### Install eckit
  install-dep.sh --owner ecmwf --repo eckit --branch master --prefix ${PREFIX} ${cleanup} --cmake "${CMAKE_OPTIONS} -DENABLE_TESTS=OFF -DENABLE_ECKIT_SQL=OFF"

  ### Install fckit
  install-dep.sh --owner ecmwf --repo fckit --branch master --prefix ${PREFIX} ${cleanup} --cmake "${CMAKE_OPTIONS} -DENABLE_TESTS=OFF"

  ### Install gridtools
  if ${with_gridtools}; then
    install-dep.sh --owner gridtools --repo gridtools --version 80187f11 --prefix ${PREFIX} ${cleanup} --cmake "${CMAKE_OPTIONS} -DBUILD_TESTING=OFF -DGT_ENABLE_OPENMP=OFF"
    export GridTools_DIR=${PREFIX}/lib/cmake # see GridTools issue (https://github.com/GridTools/gridtools/issues/1395)
  fi

else

  ### Load ecbuild
  if [[ -f ${PREFIX}/share/ecbuild/cmake/ecbuild.cmake ]]; then
    ECBUILD_MODULE_PATH=${PREFIX}/share/ecbuild/cmake
  fi

  ### Fix non-standard GridTools installation detection
  if [[ -f ${PREFIX}/lib/cmake/GridToolsConfig.cmake ]]; then
    export GridTools_DIR=${PREFIX}/lib/cmake # see GridTools issue (https://github.com/GridTools/gridtools/issues/1395)
  fi

fi

### Install atlas
SOURCE_DIR=$SCRIPTDIR/..
BUILD_DIR=${TMPDIR}/builds/atlas
mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake -DCMAKE_MODULE_PATH=${ECBUILD_MODULE_PATH} -DCMAKE_INSTALL_PREFIX=${PREFIX} ${CMAKE_OPTIONS} ${SOURCE_DIR}
make -j8
make install

#if cleanup is empty string
echo "'--cleanup' option was not specified. To clean up manually, remove content in:"
echo "    ${TMPDIR}/downloads"
echo "    ${TMPDIR}/builds"
#fi