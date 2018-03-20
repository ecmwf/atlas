#!/bin/bash

set -e

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MPI="$1"
os=$(uname)
OMPIVER=openmpi-3.0.0
MPICHVER=3.2.1

PREFIX=$(pwd)/${MPI}
mkdir -p ${PREFIX}
touch ${PREFIX}/env.sh

MPI_INSTALLED=false

case "$os" in
    Darwin)
        case "$MPI" in
            mpich)
                brew upgrade mpich || brew install mpich
                ;;
            openmpi)
                brew upgrade openmpi || brew install openmpi
                echo "localhost slots=12" >> /usr/local/etc/openmpi-default-hostfile
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
        case "$MPI" in
            mpich)
                if [ -f mpich/include/mpi.h ]; then
                  echo "mpich/include/mpi.h found."
                fi
                if [ -f mpich/lib/libmpich.so ]; then
                  echo "libmpich.so found -- nothing to build."
                else
                  echo "Downloading mpich source..."
                  wget http://www.mpich.org/static/downloads/${MPICHVER}/mpich-${MPICHVER}.tar.gz
                  tar xfz mpich-${MPICHVER}.tar.gz
                  rm mpich-${MPICHVER}.tar.gz
                  echo "Configuring and building mpich..."
                  cd mpich-${MPICHVER}
                  ${SCRIPTDIR}/reduce-output.sh ./configure \
                          --prefix=${PREFIX} \
                          --enable-static=false \
                          --enable-alloca=true \
                          --enable-threads=single \
                          --enable-fortran=yes \
                          --enable-fast=all \
                          --enable-g=none \
                          --enable-timing=none
                  ${SCRIPTDIR}/reduce-output.sh make -j4
                  ${SCRIPTDIR}/reduce-output.sh make install
                  MPI_INSTALLED=true
                  cd -
                fi
                ;;
            openmpi)
                if [ -f openmpi/include/mpi.h ]; then
                  echo "openmpi/include/mpi.h found."
                fi
                if [ -f openmpi/lib/libmpi.so ] || [ -f openmpi/lib64/libmpi.so ]; then
                  echo "libmpi.so found -- nothing to build."
                else
                  echo "Downloading openmpi source..."
                  wget --no-check-certificate https://www.open-mpi.org/software/ompi/v3.0/downloads/$OMPIVER.tar.gz
                  tar -zxf $OMPIVER.tar.gz
                  rm $OMPIVER.tar.gz
                  echo "Configuring and building openmpi..."
                  cd $OMPIVER
                  ${SCRIPTDIR}/reduce-output.sh ./configure \
                          --prefix=${PREFIX}
                  ${SCRIPTDIR}/reduce-output.sh make -j4
                  ${SCRIPTDIR}/reduce-output.sh make install
                  MPI_INSTALLED=true
                  cd -
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
fi

echo "Please source ${PREFIX}/env.sh, containing:"
cat ${PREFIX}/env.sh
