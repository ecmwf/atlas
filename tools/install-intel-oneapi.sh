#!/bin/sh

KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
wget https://apt.repos.intel.com/intel-gpg-keys/$KEY
sudo apt-key add $KEY
rm $KEY
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update
sudo apt-get install \
    intel-oneapi-compiler-fortran \
    intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic \
    intel-oneapi-mpi \
    intel-oneapi-mpi-devel \
    intel-oneapi-mkl
