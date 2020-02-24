#!/usr/bin/env bash
HERE="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
source ${HERE}/source-me.sh
build hello_world_fortran project_hello_world_fortran
