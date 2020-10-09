#! /usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -e -o pipefail
set +x

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

# Some defaults for the arguments
owner=ecmwf
repo="unknown"
branch=master
PREFIX=$(pwd)/install
skip_install=true
cleanup=false
version=${branch}

while [ $# != 0 ]; do
    case "$1" in
    "--owner")
        owner="$2"; shift
        ;;
    "--repo")
        repo="$2"; shift
        ;;
    "--branch")
        version="$2"; shift
        ;;
    "--version")
        version="$2"; shift
        ;;
    "--depends")
        depends=($2); shift
        ;;
    "--cmake")
        CMAKE_OPTIONS="$2"; shift
        ;;
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--no-cache")
        skip_install=false;
        ;;
    "--cleanup")
        cleanup=true;
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

function cleanup_atexit {
  EXIT_CODE=$?
  { set +ex; } 2>/dev/null
  if [ $EXIT_CODE != 0 ]; then
    set -x
    rm -f ${PREFIX}/${repo}-sha-${SHA}
    rm -f ${PREFIX}/sha/${repo}
  fi
  exit $EXIT_CODE
}
trap cleanup_atexit EXIT

SHA=$(github-sha.sh ${owner} ${repo} ${version})
echo "Installing ${owner}/${repo} version ${version} at sha [${SHA}]"

all_deps_uptodate=true
for dep in "${depends[@]}"; do
  DEP=$(echo ${dep} | awk '{print toupper($0)}')
  DEP_PATH=$(eval "echo \$${DEP}_PATH")
  curr_sha=$(cat ${DEP_PATH}/sha/${dep})

  if [[ -f ${PREFIX}/sha/${dep} ]]; then
    prev_sha=$(cat ${PREFIX}/sha/${dep})
    if [ "${prev_sha}" == "${curr_sha}" ]; then
      uptodate=true
    else
      uptodate=false
    fi
  else
    uptodate=true
  fi
  mkdir -p ${PREFIX}/sha
  touch ${PREFIX}/sha/${dep}
  echo ${curr_sha} > ${PREFIX}/sha/${dep}

  if ! ${uptodate}; then
    all_deps_uptodate=false
  fi
done
if ${skip_install} ; then
  skip_install=${all_deps_uptodate}
fi

if [[ -f "${PREFIX}/${repo}-sha-${SHA}" ]]; then
  echo "Cached ${repo} (${repo}-sha-${SHA}) is up to date."
else
  skip_install=false
  if ls ${PREFIX}/${repo}-sha-* 1> /dev/null 2>&1; then
     echo "Cached ${repo} $(ls -1 ${PREFIX}/${repo}-sha-*) is out of date."
  fi
fi

if ! ${skip_install} ; then
  echo "Installing ${repo} ${repo}-sha-${SHA}"
  pushd . >> /dev/null 2>&1
  if [ -z "${TMPDIR}" ]; then
    TMPDIR=${HOME}/tmp
  fi
  SOURCE_DIR=${TMPDIR}/downloads/${repo}
  BUILD_DIR=${TMPDIR}/builds/${repo}
  rm -rf ${SOURCE_DIR} ${BUILD_DIR}

  # Download
  mkdir -p ${TMPDIR}/downloads
  echo "+ git clone https://github.com/${owner}/${repo} ${SOURCE_DIR}"
  git clone https://github.com/${owner}/${repo} ${SOURCE_DIR}
  echo "+ cd ${SOURCE_DIR}"
  cd ${SOURCE_DIR}
  echo "+ git checkout ${version}"
  git checkout ${version}
  
  # Configure and install
  mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
  if [ -z ${ECBUILD_MODULE_PATH+x} ]; then
    echo "No ECBUILD_MODULE_PATH defined"
  else
    SET_CMAKE_MODULE_PATH="-DCMAKE_MODULE_PATH=${ECBUILD_MODULE_PATH}"
  fi
  echo "+ cmake ${SET_CMAKE_MODULE_PATH} -DCMAKE_INSTALL_PREFIX=${PREFIX} ${CMAKE_OPTIONS} ${SOURCE_DIR}"
  cmake ${SET_CMAKE_MODULE_PATH} -DCMAKE_INSTALL_PREFIX=${PREFIX} ${CMAKE_OPTIONS} ${SOURCE_DIR}
  core_count=$(nproc || echo 4)
  echo "+ make -j${core_count} install"
  make -j${core_count} install

  # Store version information
  mkdir -p ${PREFIX}/sha
  touch ${PREFIX}/${repo}-sha-${SHA}
  echo ${SHA} > ${PREFIX}/sha/${repo}
  for dep in "${depends[@]}"; do
    DEP=$(echo ${dep} | awk '{print toupper($0)}')
    DEP_PATH=$(eval "echo \$${DEP}_PATH")
    curr_sha=$(cat ${DEP_PATH}/sha/${dep})
    echo ${curr_sha} > ${PREFIX}/sha/${dep}
  done
  popd >> /dev/null 2>&1
  if ${cleanup} ; then
    rm -rf ${SOURCE_DIR} ${BUILD_DIR}
  fi
fi

# Export
#REPO=$(echo $repo | awk '{print toupper($0)}')
#echo "To make it easier to find ${repo}, execute now"
#echo "export ${REPO}_PATH=${PREFIX}"
