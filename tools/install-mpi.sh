#!/bin/bash


set +x
set -e -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

# Some defaults for the arguments
PREFIX=$(pwd)/${MPI}
mpi_override=false
MPI=openmpi

while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        PREFIX="$2"; shift
        ;;
    "--override")
        mpi_override=true; 
        ;;
    "--version")
        mpi_version="$2"; shift
        ;;
    "--mpi")
        MPI="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

os=$(uname)
OMPIVER=4.1.1
MPICHVER=3.4.2

if [ ! -z ${mpi_version+x} ]; then
  if [[ "${MPI}" =~ [Oo][Pp][Ee][Nn]\-?[Mm][Pp][Ii] ]]; then
    OMPIVER=${mpi_version}
  fi
  if [[ "${MPI}" =~ [Mm][Pp][Ii][Cc][Hh] ]]; then
    MPICHVER=${mpi_version}
  fi
fi


mkdir -p ${PREFIX}
touch ${PREFIX}/env.sh

MPI_INSTALLED=false

case "$os" in
    Darwin)
        case "$MPI" in
            mpich)
                brew ls --versions mpich || brew install mpich
                ;;
            openmpi)
                brew ls --versions openmpi || brew install openmpi
                echo "localhost slots=72" >> $(brew --prefix)/etc/openmpi-default-hostfile
                # workaround for open-mpi/omp#7516
                echo "setting the mca gds to hash..."
                echo "gds = hash" >> $(brew --prefix)/etc/pmix-mca-params.conf

                # workaround for open-mpi/ompi#5798
                echo "setting the mca btl_vader_backing_directory to /tmp..."
                echo "btl_vader_backing_directory = /tmp" >> $(brew --prefix)/etc/openmpi-mca-params.conf
                ;;
            *)
                echo "Unknown MPI implementation: $MPI"
                exit 1
                ;;
        esac
    ;;

    Linux)
        if [ -n "${MPI_HOME}" ]; then
          echo "MPI is already installed at MPI_HOME=${MPI_HOME}."
          echo "Not taking any action."
          exit 0
        fi
        if [ -n "${I_MPI_ROOT}" ]; then
          echo "MPI is already installed at I_MPI_ROOT=${I_MPI_ROOT}."
          echo "Not taking any action."
          exit 0
        fi
        case "$MPI" in
            mpich)
                if [ -f ${PREFIX}/include/mpi.h ]; then
                  echo "${PREFIX}/include/mpi.h found"
                fi
                if [ -f ${PREFIX}/lib/libmpich.so ]; then
                  echo "${PREFIX}/lib/libmpich.so found -- nothing to build."
                else
                  echo "Downloading mpich source..."
                  wget http://www.mpich.org/static/downloads/${MPICHVER}/mpich-${MPICHVER}.tar.gz
                  tar xfz mpich-${MPICHVER}.tar.gz
                  rm mpich-${MPICHVER}.tar.gz
                  echo "Configuring and building mpich..."
                  cd mpich-${MPICHVER}
                  unset F90
                  unset F90FLAGS
                  ${SCRIPTDIR}/reduce-output.sh ./configure \
                          --prefix=${PREFIX} \
                          --enable-static=false \
                          --enable-alloca=true \
                          --enable-threads=single \
                          --enable-fortran=yes \
                          --enable-fast=all \
                          --enable-g=none \
                          --enable-timing=none
                  ${SCRIPTDIR}/reduce-output.sh make -j48
                  ${SCRIPTDIR}/reduce-output.sh make install
                  MPI_INSTALLED=true
                  cd -
                  rm -rf mpich-${MPICHVER}
                fi
                ;;
            openmpi)
                if [ -f ${PREFIX}/include/mpi.h ]; then
                  echo "openmpi/include/mpi.h found."
                fi
                if [ -f ${PREFIX}/lib/libmpi.so ] || [ -f ${PREFIX}/lib64/libmpi.so ]; then
                  echo "libmpi.so found -- nothing to build."
                else
                  echo "Downloading openmpi source..."
                  wget --no-check-certificate https://www.open-mpi.org/software/ompi/v4.1/downloads/openmpi-$OMPIVER.tar.gz
                  tar -zxf openmpi-$OMPIVER.tar.gz
                  rm openmpi-$OMPIVER.tar.gz
                  echo "Configuring and building openmpi..."
                  cd openmpi-$OMPIVER
                  ${SCRIPTDIR}/reduce-output.sh ./configure --prefix=${PREFIX}
                  ${SCRIPTDIR}/reduce-output.sh make -j4
                  ${SCRIPTDIR}/reduce-output.sh make install
                  MPI_INSTALLED=true
                  echo "localhost slots=72" >> ${PREFIX}/etc/openmpi-default-hostfile
                  cd -
                  rm -rf openmpi-$OMPIVER
                fi
                ;;
            *)
                echo "Unknown MPI implementation: $MPI"
                exit 1
                ;;
        esac
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac


if ${MPI_INSTALLED} ; then
cat > ${PREFIX}/env.sh << EOF
export MPI_HOME=${PREFIX}
export PATH=\${MPI_HOME}/bin:\${PATH}
EOF
echo "Please source ${PREFIX}/env.sh, containing:"
cat ${PREFIX}/env.sh
fi
