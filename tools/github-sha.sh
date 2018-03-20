#!/usr/bin/env bash

owner=$1
repo=$2
branch=$3

if [ -z "${TMPDIR}" ]; then
  TMPDIR=${HOME}/tmp
fi
SOURCE_DIR=${TMPDIR}/check-git-${repo}

dump_output() {
   echo " ++ Tailing the last 500 lines of output from ${BUILD_OUTPUT}"
   tail -500 ${BUILD_OUTPUT}
}
error_handler() {
  echo ERROR: An error was encountered with the build.
  rm -rf ${SOURCE_DIR}
  dump_output
  exit 1
}
# If an error occurs, run our error handler to output a tail of the build
trap 'error_handler' ERR

BUILD_OUTPUT=${TMPDIR}/tmp-${repo}.log

mkdir -p ${TMPDIR}
touch $BUILD_OUTPUT
pushd . >> ${BUILD_OUTPUT} 2>&1
git clone --depth=1 -b ${branch} https://github.com/${owner}/${repo} ${SOURCE_DIR} >> ${BUILD_OUTPUT} 2>&1 
cd ${SOURCE_DIR} >> ${BUILD_OUTPUT} 2>&1 
git rev-parse HEAD
cd ..
rm -rf ${SOURCE_DIR}
popd >> ${BUILD_OUTPUT} 2>&1
