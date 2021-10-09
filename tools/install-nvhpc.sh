#!/bin/sh

# Install NVHPC
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

version=21.9

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export NVHPC_INSTALL_DIR=$(pwd)/nvhpc-install
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
    "--version")
        version="$2"; shift
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

if [ -d "${NVHPC_INSTALL_DIR}" ]; then
    if [[ $(find "${NVHPC_INSTALL_DIR}" -name "nvc" | wc -l) == 1 ]]; then
      echo "NVHPC already installed at ${NVHPC_INSTALL_DIR}"
      exit
    fi
fi

# Example download URL for version 21.9
#    https://developer.download.nvidia.com/hpc-sdk/21.9/nvhpc_2020_219_Linux_x86_64_cuda_11.0.tar.gz

ver="$(echo $version | tr -d . )"
URL=$(curl -s "https://developer.nvidia.com/nvidia-hpc-sdk-$ver-downloads" | grep -oP "https://developer.download.nvidia.com/hpc-sdk/([0-9]{2}\.[0-9]+)/nvhpc_([0-9]{4})_([0-9]+)_Linux_$(uname -m)_cuda_([0-9\.]+).tar.gz" | sort | tail -1)
FOLDER="$(basename "$(echo "${URL}" | grep -oP '[^/]+$')" .tar.gz)"

if [ ! -d "${TEMPORARY_FILES}/${FOLDER}" ]; then
  echo "Downloading ${TEMPORARY_FILES}/${FOLDER} from URL [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  curl --location \
       --user-agent "pgi-travis (https://github.com/nemequ/pgi-travis)" \
       "${URL}" | tar zx -C "${TEMPORARY_FILES}"
else
   echo "Download already present in ${TEMPORARY_FILES}/${FOLDER}"
fi

echo "+ ${TEMPORARY_FILES}/${FOLDER}/install"
"${TEMPORARY_FILES}/${FOLDER}/install"

#comment out to cleanup
#rm -rf "${TEMPORARY_FILES}/${FOLDER}"

NVHPC_VERSION=$(basename "${NVHPC_INSTALL_DIR}"/Linux_$(uname -m)/*.*/)

# Use gcc which is available in PATH
${NVHPC_INSTALL_DIR}/Linux_$(uname -m)/${NVHPC_VERSION}/compilers/bin/makelocalrc \
  -x ${NVHPC_INSTALL_DIR}/Linux_$(uname -m)/${NVHPC_VERSION}/compilers/bin \
  -gcc $(which gcc) \
  -gpp $(which g++) \
  -g77 $(which gfortran)

cat > ${NVHPC_INSTALL_DIR}/env.sh << EOF
### Variables
export NVHPC_INSTALL_DIR=${NVHPC_INSTALL_DIR}
export NVHPC_VERSION=${NVHPC_VERSION}
export NVHPC_DIR=\${NVHPC_INSTALL_DIR}/Linux_$(uname -m)/\${NVHPC_VERSION}

### Compilers
export PATH=\${NVHPC_DIR}/compilers/bin:\${PATH}
export NVHPC_LIBRARY_PATH=\${NVHPC_DIR}/compilers/lib
export LD_LIBRARY_PATH=\${NVHPC_LIBRARY_PATH}

### MPI
export MPI_HOME=\${NVHPC_DIR}/comm_libs/mpi
export PATH=\${MPI_HOME}/bin:\${PATH}
EOF

cat ${NVHPC_INSTALL_DIR}/env.sh

