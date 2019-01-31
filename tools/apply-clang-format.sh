#!/usr/bin/env bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPTDIR/../src

echo "Applying $(clang-format --version) ..."

find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file

echo "Applying $(clang-format --version) ... done"

