#!/usr/bin/env bash

# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


set -e

function help {
echo "$(basename ${BASH_SOURCE[0]}) --project PROJECT_DIR [--cxx] [--fortran] [--build]"
echo ""
echo "DECRIPTION"
echo "  Tool that generates a simple CMake project (not with ecbuild)"
echo "  with a \"main\" program that is linked to atlas."
echo "  Note: It is intended for quickly prototyping or testing atlas features."
echo "        For more complicated projects, please see"
echo "        https://sites.ecmwf.int/docs/atlas/getting_started/linking"
echo ""
echo "OPTIONS"
echo "  --help              Print this help"
echo "  --project           Path of project directory to be generated"
echo "  --c++ , --cxx       Create C++ project (default if not specified)"
echo "  --fortran           Create Fortran project instead of C++"
echo "  --build             Attempt to build template project in 'build'"
echo "                      subdirectory of generated project"
}

LANGUAGES=C
WITH_CXX=false
WITH_Fortran=false
WITH_build=false

if [ $# -eq 0 ]; then
    help
    exit 1
fi

while test $# -gt 0; do

  case "$1" in
    --project)
    PROJECT_DIR=$2
    PROJECT_NAME=$(basename $PROJECT_DIR)
    shift # past argument
    shift # past value 
    ;;
    --c++|--cxx)
    WITH_CXX=true
    LANGUAGES="${LANGUAGES} CXX"
    shift
    ;;
    --fortran)
    WITH_Fortran=true
    LANGUAGES="${LANGUAGES} Fortran"
    shift
    ;;
    --build)
    WITH_build=true
    shift
    ;;
    --help|-h)
    help
    shift
    ;;
    *)
    # unknown option
    help
    break
    ;;
  esac
done

if ( $WITH_Fortran && $WITH_CXX ) ; then
  echo "ERROR: Both --cxx and --fortran were supplied. Choose one only."
  exit 1
fi

if ( ! $WITH_Fortran && ! $WITH_CXX ) ; then
  echo "No language specified. Defaulting to use C++."
  WITH_CXX=true
  LANGUAGES="${LANGUAGES} CXX"
fi

mkdir -p $PROJECT_DIR
cd $PROJECT_DIR
cat << EOF > CMakeLists.txt
cmake_minimum_required( VERSION 3.18 )

project( ${PROJECT_NAME} VERSION 0.0.0 LANGUAGES ${LANGUAGES} )

find_package( atlas REQUIRED )

EOF

if ( $WITH_CXX ) ; then
  cat <<- EOF >> CMakeLists.txt
	add_executable( main )
	target_sources( main PUBLIC main.cc )
	target_link_libraries( main PUBLIC atlas )
	EOF
  cat <<- EOF > main.cc 
	#include "atlas/library.h"
	#include "atlas/runtime/Log.h"
	int main(int argc, char* argv[]) {
	    atlas::initialize(argc,argv);
	    atlas::Log::info() << "Hello from C++" << std::endl; 
	    atlas::finalize();
	    return 0;
	}
	EOF
elif ( $WITH_Fortran ) ; then
  cat <<- EOF >> CMakeLists.txt
	add_executable( main )
	target_sources( main PUBLIC main.F90 )
	target_link_libraries( main PUBLIC atlas_f )
	EOF
  cat <<- EOF > main.F90
	program main
	use atlas_module
	implicit none
	call atlas_initialize()
	call atlas_log%info("Hello from Fortran")
	call atlas_finalize()
	end program
	EOF
fi

echo "Template project ${PROJECT_NAME} has been created in ${PROJECT_DIR}"

if ( $WITH_build ) ; then
    mkdir -p build
    cd build

    echo ""
    echo "Project ${PROJECT_NAME} is being configured in $(pwd)"

    cmake ..

    echo ""
    echo "Project ${PROJECT_NAME} is being built in $(pwd)"

    make VERBOSE=1
    
    echo ""
    echo "Running program main in $(pwd)"
    echo ""
    ./main
fi
