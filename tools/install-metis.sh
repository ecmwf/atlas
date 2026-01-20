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
metis_version=5.2.1
metis_version=master
gklib_version=master

while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--version")
        metis_version="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done


echo "Installing metis version ${metis_version}"

metis_installed=${PREFIX}/metis-${metis_version}-installed
if [[ -f "${metis_installed}" ]]; then
  echo "metis ${metis_version} is already installed at ${PREFIX}"
  exit
fi

os=$(uname)
case "$os" in
    Darwin)
      brew install metis
      exit
    ;;
    *)
    ;;
esac


if [ -z "${TMPDIR+x}" ]; then
  TMPDIR=${HOME}/tmp
fi
mkdir -p ${TMPDIR}/downloads

gklib_tag=${gklib_version}
gklib_tarball_url=https://github.com/KarypisLab/GKlib/archive/refs/heads/master.tar.gz
gklib_tarball=$TMPDIR/downloads/gklib-${gklib_tag}.tar.gz
gklib_dir=$TMPDIR/downloads/gklib-${gklib_version}

if [ ! -d "${gklib_dir}" ]; then
  echo "+ curl -L ${gklib_tarball_url} > ${gklib_tarball}"
  curl -L ${gklib_tarball_url} > ${gklib_tarball}
  echo "+ tar xzf ${gklib_tarball} -C ${TMPDIR}/downloads"
  tar xzf ${gklib_tarball} -C ${TMPDIR}/downloads
fi
echo "+ cd ${gklib_dir}"
cd ${gklib_dir}
echo "Applying patch to remove so-versioning from GKlib"
sed -i -e 's/^set_target_properties(${PROJECT_NAME} PROPERTIES/\# set_target_properties(${PROJECT_NAME} PROPERTIES/g' CMakeLists.txt
sed -i -e 's/^  SOVERSION ${PROJECT_VERSION_MAJOR}/\#   SOVERSION ${PROJECT_VERSION_MAJOR}/g' CMakeLists.txt
sed -i -e 's/^  VERSION   ${PROJECT_VERSION})/\#   VERSION   ${PROJECT_VERSION})/g' CMakeLists.txt
make config shared=1 prefix=${PREFIX}
make -j8
make install

metis_tag=v${metis_version}
metis_tarball_url=https://github.com/KarypisLab/METIS/archive/refs/tags/v${metis_version}.tar.gz

if [[ "${metis_version}" == "master" ]]; then
  metis_tag=${metis_version}
  metis_tarball_url=https://github.com/KarypisLab/METIS/archive/refs/heads/master.tar.gz
fi

metis_tarball=$TMPDIR/downloads/metis-${metis_tag}.tar.gz
metis_dir=$TMPDIR/downloads/metis-${metis_version}

if [ ! -d "${metis_dir}" ]; then
  echo "+ curl -L ${metis_tarball_url} > ${metis_tarball}"
  curl -L ${metis_tarball_url} > ${metis_tarball}
  echo "+ tar xzf ${metis_tarball} -C ${TMPDIR}/downloads"
  tar xzf ${metis_tarball} -C ${TMPDIR}/downloads
fi
echo "+ cd ${metis_dir}"
cd ${metis_dir}
echo "Applying patch to link metis to GKlib"
sed -i -e 's/add_library(metis ${METIS_LIBRARY_TYPE} ${metis_sources})/add_library(metis ${METIS_LIBRARY_TYPE} ${metis_sources})\ntarget_link_libraries(metis PUBLIC GKlib)/g' libmetis/CMakeLists.txt
make config shared=1 prefix=${PREFIX} gklib_path=${PREFIX}
echo "+ make -j8"
make -j8
echo "+ make install"
make install

echo "+ rm -rf \${gklib_tarball} \${gklib_dir}"
rm -rf ${gklib_tarball} ${gklib_dir}
echo "+ rm -rf \${metis_tarball} \${metis_dir}"
rm -rf ${metis_tarball} ${metis_dir}

echo "+ touch ${metis_installed}"
touch ${metis_installed}
