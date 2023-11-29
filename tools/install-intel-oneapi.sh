#!/usr/bin/env bash

version=2023.2.0
KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
wget https://apt.repos.intel.com/intel-gpg-keys/$KEY
sudo apt-key add $KEY
rm $KEY
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update
sudo apt-get install \
    intel-oneapi-compiler-fortran-$version \
    intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-$version \
    intel-oneapi-mpi-devel-2021.10.0 \
    intel-oneapi-mkl-$version
