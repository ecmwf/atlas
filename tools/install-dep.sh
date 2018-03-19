#! /usr/bin/env bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

# Some defaults for the arguments
owner=ecmwf
repo="unknown"
branch=master
PREFIX=$(pwd)/install
skip_install=true

while [ $# != 0 ]; do
    case "$1" in
    "--owner")
        owner="$2"; shift
        ;;
    "--repo")
        repo="$2"; shift
        ;;
    "--branch")
        branch="$2"; shift
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
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

SHA=$(github-sha.sh ${owner} ${repo} ${branch})
echo "Installing ${owner}/${repo} branch ${branch} at sha [${SHA}]"

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

if [[ -f "${PREFIX}/sha-${SHA}" ]]; then
  echo "Cached ${repo} (sha-${SHA}) is up to date."
else
  skip_install=false
  if ls ${PREFIX}/sha-* 1> /dev/null 2>&1; then
     echo "Cached ${repo} $(ls -1 ${PREFIX}/sha-*) is out of date."
  fi
fi

if ! ${skip_install} ; then
  echo "Installing ${repo} sha-${SHA}"
  pushd . >> /dev/null 2>&1
  if [ -z "${TMPDIR}" ]; then
    TMPDIR=${HOME}/tmp
  fi
  SOURCE_DIR=${TMPDIR}/downloads/${repo}
  BUILD_DIR=${TMPDIR}/builds/${repo}
  rm -rf ${PREFIX} ${SOURCE_DIR} ${BUILD_DIR}

  # Download
  mkdir -p ${TMPDIR}/downloads
  echo "+ git clone -b ${branch} https://github.com/${owner}/${repo} ${SOURCE_DIR}"
  git clone -b ${branch} https://github.com/${owner}/${repo} ${SOURCE_DIR}

  # Configure and install
  mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
  echo "+ cmake -DCMAKE_MODULE_PATH=${ECBUILD_MODULE_PATH} -DCMAKE_INSTALL_PREFIX=${PREFIX} ${CMAKE_OPTIONS} ${SOURCE_DIR}"
  cmake -DCMAKE_MODULE_PATH=${ECBUILD_MODULE_PATH} -DCMAKE_INSTALL_PREFIX=${PREFIX} ${CMAKE_OPTIONS} ${SOURCE_DIR}
  make -j4 install

  # Store version information
  mkdir -p ${PREFIX}/sha
  touch ${PREFIX}/sha-${SHA}
  echo ${SHA} > ${PREFIX}/sha/${repo}
  for dep in "${depends[@]}"; do
    DEP=$(echo ${dep} | awk '{print toupper($0)}')
    DEP_PATH=$(eval "echo \$${DEP}_PATH")
    curr_sha=$(cat ${DEP_PATH}/sha/${dep})
    echo ${curr_sha} > ${PREFIX}/sha/${dep}
  done
  popd >> /dev/null 2>&1
  rm -rf ${SOURCE_DIR} ${BUILD_DIR}

fi

# Export
REPO=$(echo $repo | awk '{print toupper($0)}')
echo "To make it easier to find ${repo}, execute now"
echo "export ${REPO}_PATH=${PREFIX}"
