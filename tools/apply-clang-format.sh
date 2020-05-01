#!/usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

_REQUIRED_CLANG_VERSION='9.0.1'

export PATH=/usr/local/apps/clang/9.0.1/bin:$PATH

if ! [ -x "$(command -v clang-format)" ]; then
  echo 'Error: clang-format is not installed.' >&2
  exit 1
fi

if ! [[ $(clang-format --version) =~ ${_REQUIRED_CLANG_VERSION} ]]; then
    echo "Error: Require clang-format version: ${_REQUIRED_CLANG_VERSION}"
    echo "    > $(which clang-format) --version"
    echo "      $(clang-format --version)"
    if [[ $1 =~ --no-version-check ]]; then
        echo "--no-version-check --> Continue anyway"
    else
        exit 1
    fi
fi

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPTDIR/../src

echo "Applying $(clang-format --version) ..."

find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file

echo "Applying $(clang-format --version) ... done"

