#!/usr/bin/env bash
HERE="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
source ${HERE}/source-me.sh
build hello_world project_hello_world
