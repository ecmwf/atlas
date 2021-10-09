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
fftw_version=3.3.10
fftw_configure="--enable-shared"
fftw_with_single=false

while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--version")
        fftw_version="$2"; shift
        ;;
    "--with-single")
        fftw_with_single=true;
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

echo "Installing FFTW version ${fftw_version}"

fftw_installed=${PREFIX}/fftw-${fftw_version}-installed
if [[ -f "${fftw_installed}" ]]; then
  echo "FFTW ${fftw_version} is already installed at ${PREFIX}"
  exit
fi

os=$(uname)
case "$os" in
    Darwin)
      brew ls --versions fftw || brew install fftw
      exit
    ;;
    *)
    ;;
esac


if [ -z "${TMPDIR+x}" ]; then
  TMPDIR=${HOME}/tmp
fi
mkdir -p ${TMPDIR}/downloads

fftw_tarball_url=http://www.fftw.org/fftw-${fftw_version}.tar.gz
fftw_tarball=$TMPDIR/downloads/fftw-${fftw_version}.tar.gz
fftw_dir=$TMPDIR/downloads/fftw-${fftw_version}

echo "+ curl -L ${fftw_tarball_url} > ${fftw_tarball}"
curl -L ${fftw_tarball_url} > ${fftw_tarball}
echo "+ tar xzf ${fftw_tarball} -C ${TMPDIR}/downloads"
tar xzf ${fftw_tarball} -C ${TMPDIR}/downloads
echo "+ cd ${fftw_dir}"
cd ${fftw_dir}
echo "+ ./configure --prefix=${PREFIX} ${fftw_configure}"
./configure --prefix=${PREFIX} ${fftw_configure}
echo "+ make -j8"
make -j8
echo "+ make install"
make install

if $fftw_with_single; then
  # Now again in single precision
  make clean
  echo "+ ./configure --prefix=${PREFIX} ${fftw_configure} --enable-float"
  ./configure --prefix=${PREFIX} ${fftw_configure} --enable-float
  echo "+ make -j8"
  make -j8
  echo "+ make install"
  make install
fi


echo "+ rm -rf \${fftw_tarball} \${fftw_dir}"
rm -rf ${fftw_tarball} ${fftw_dir}

echo "+ touch ${fftw_installed}"
touch ${fftw_installed}
