#!/usr/bin/env bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPTDIR/../src
find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file

