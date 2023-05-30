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
qhull_version=8.1-alpha3

while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--version")
        qhull_version="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

echo "Installing Qhull version ${qhull_version}"

qhull_installed=${PREFIX}/qhull-${qhull_version}-installed
if [[ -f "${qhull_installed}" ]]; then
  echo "Qhull ${qhull_version} is already installed at ${PREFIX}"
  exit
fi

os=$(uname)
case "$os" in
    Darwin)
      brew install qhull
      exit
    ;;
    *)
    ;;
esac


if [ -z "${TMPDIR+x}" ]; then
  TMPDIR=${HOME}/tmp
fi
mkdir -p ${TMPDIR}/downloads

qhull_tarball_url=https://github.com/qhull/qhull/archive/refs/tags/v${qhull_version}.tar.gz
qhull_tarball=$TMPDIR/downloads/qhull-v${qhull_version}.tar.gz
qhull_dir=$TMPDIR/downloads/qhull-${qhull_version}

echo "+ curl -L ${qhull_tarball_url} > ${qhull_tarball}"
curl -L ${qhull_tarball_url} > ${qhull_tarball}
echo "+ tar xzf ${qhull_tarball} -C ${TMPDIR}/downloads"
tar xzf ${qhull_tarball} -C ${TMPDIR}/downloads
echo "+ cd ${qhull_dir}"
cd ${qhull_dir}
#echo "+ ./configure --prefix=${PREFIX} ${fftw_configure}"
#./configure --prefix=${PREFIX} ${fftw_configure}
echo "+ mkdir build-release"
mkdir build-release
echo "+ cd build-release"
cd build-release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX}
echo "+ make -j8"
make -j8
echo "+ make install"
make install

echo "+ rm -rf \${qhull_tarball} \${qhull_dir}"
rm -rf ${qhull_tarball} ${qhull_dir}

echo "+ touch ${qhull_installed}"
touch ${qhull_installed}
