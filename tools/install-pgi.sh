#!/bin/sh

# Install PGI Community Edition on Travis
# https://github.com/nemequ/pgi-travis
#
# Originally written for Squash <https://github.com/quixdb/squash> by
# Evan Nemerson.  For documentation, bug reports, support requests,
# etc. please use <https://github.com/nemequ/pgi-travis>.
#
# To the extent possible under law, the author(s) of this script have
# waived all copyright and related or neighboring rights to this work.
# See <https://creativecommons.org/publicdomain/zero/1.0/> for
# details.

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export NVHPC_INSTALL_DIR=$(pwd)/pgi-install
export NVHPC_SILENT=true
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export NVHPC_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--verbose")
        export NVHPC_SILENT=false;
        ;;
    *)   
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac 
    shift
done

case "$(uname -m)" in
	x86_64|ppc64le|aarch64)
		;;
	*)
		echo "Unknown architecture: $(uname -m)" >&2
		exit 1
		;;
esac

URL="$(curl -s 'https://developer.nvidia.com/nvidia-hpc-sdk-download' | grep -oP "https://developer.download.nvidia.com/hpc-sdk/([0-9]{2}\.[0-9]+)/nvhpc_([0-9]{4})_([0-9]+)_Linux_$(uname -m)_cuda_([0-9\.]+).tar.gz" | sort | tail -1)"
FOLDER="$(basename "$(echo "${URL}" | grep -oP '[^/]+$')" .tar.gz)"

if [ ! -z "${TRAVIS_REPO_SLUG}" ]; then
	curl --location \
    	--user-agent "pgi-travis (https://github.com/nemequ/pgi-travis; ${TRAVIS_REPO_SLUG})" \
    	--referer "http://www.pgroup.com/products/community.htm" \
    	--header "X-Travis-Build-Number: ${TRAVIS_BUILD_NUMBER}" \
    	--header "X-Travis-Event-Type: ${TRAVIS_EVENT_TYPE}" \
    	--header "X-Travis-Job-Number: ${TRAVIS_JOB_NUMBER}" \
    	"${URL}" | tar zx -C "${TEMPORARY_FILES}"
else

  if [ ! -d "${TEMPORARY_FILES}/${FOLDER}" ]; then
    echo "Downloading ${TEMPORARY_FILES}/${FOLDER} from URL [${URL}]"
    mkdir -p ${TEMPORARY_FILES}
    curl --location \
         --user-agent "pgi-travis (https://github.com/nemequ/pgi-travis)" \
         "${URL}" | tar zx -C "${TEMPORARY_FILES}"
  else
     echo "Download already present in ${TEMPORARY_FILES}/${FOLDER}"
  fi

fi

echo "+ ${TEMPORARY_FILES}/${FOLDER}/install"
"${TEMPORARY_FILES}/${FOLDER}/install"

#comment out to cleanup
#rm -rf "${TEMPORARY_FILES}/${FOLDER}"

PGI_VERSION=$(basename "${NVHPC_INSTALL_DIR}"/Linux_$(uname -m)/*.*/)

# Use gcc which is available in PATH
${NVHPC_INSTALL_DIR}/Linux_$(uname -m)/${PGI_VERSION}/compilers/bin/makelocalrc \
  -x ${NVHPC_INSTALL_DIR}/Linux_$(uname -m)/${PGI_VERSION}/compilers/bin \
  -gcc $(which gcc) \
  -gpp $(which g++) \
  -g77 $(which gfortran)

cat > ${NVHPC_INSTALL_DIR}/env.sh << EOF
### Variables
PGI_INSTALL_DIR=${NVHPC_INSTALL_DIR}
PGI_VERSION=${PGI_VERSION}

### Compilers
PGI_DIR=\${PGI_INSTALL_DIR}/Linux_$(uname -m)/\${PGI_VERSION}
export PATH=\${PGI_DIR}/compilers/bin:\${PATH}
EOF

cat >> ${NVHPC_INSTALL_DIR}/env.sh << EOF

### MPI
export MPI_HOME=\${PGI_DIR}/comm_libs/mpi
export PATH=\${MPI_HOME}/bin:\${PATH}
export LD_LIBRARY_PATH=\${PGI_DIR}/compilers/lib:\${LD_LIBRARY_PATH}
EOF

