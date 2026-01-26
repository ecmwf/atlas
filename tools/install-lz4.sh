#! /usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set +x
set -e -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

# Some defaults for the arguments
PREFIX=$(pwd)/install
lz4_version=1.10.0

while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--version")
        lz4_version="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

echo "Installing lz4 version ${lz4_version}"

lz4_installed=${PREFIX}/lz4-${lz4_version}-installed
if [[ -f "${lz4_installed}" ]]; then
  echo "lz4 ${lz4_version} is already installed at ${PREFIX}"
  exit
fi

os=$(uname)
case "$os" in
    Darwin)
      brew install lz4
      exit
    ;;
    *)
    ;;
esac


if [ -z "${TMPDIR+x}" ]; then
  TMPDIR=${HOME}/tmp
fi
mkdir -p ${TMPDIR}/downloads

lz4_tarball_url=https://github.com/lz4/lz4/archive/refs/tags/v${lz4_version}.tar.gz
lz4_tarball=$TMPDIR/downloads/lz4-v${lz4_version}.tar.gz
lz4_dir=$TMPDIR/downloads/lz4-${lz4_version}

echo "+ curl -L ${lz4_tarball_url} > ${lz4_tarball}"
curl -L ${lz4_tarball_url} > ${lz4_tarball}
echo "+ tar xzf ${lz4_tarball} -C ${TMPDIR}/downloads"
tar xzf ${lz4_tarball} -C ${TMPDIR}/downloads
echo "+ cd ${lz4_dir}"
cd ${lz4_dir}
#echo "+ ./configure --prefix=${PREFIX} ${fftw_configure}"
#./configure --prefix=${PREFIX} ${fftw_configure}
echo "+ mkdir build-release"
mkdir build-release
echo "+ cd build-release"
cd build-release
cmake ../build/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX}
echo "+ make -j8"
make -j8
echo "+ make install"
make install

echo "+ rm -rf \${lz4_tarball} \${lz4_dir}"
rm -rf ${lz4_tarball} ${lz4_dir}

echo "+ touch ${lz4_installed}"
touch ${lz4_installed}
