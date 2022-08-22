#!/usr/bin/env bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

_REQUIRED_CLANG_VERSION='13.0.0'

function help() {
    cat <<EOF
USAGE:

    apply-clang-format.sh [--help] [--force] [--dryrun] [sources]

DESCRIPTION:

    Apply clang-format to source files given the atlas style.
    If no sources are given, a search for all files within src/ will be
    performed with extensions .h and .cc
EOF
  exit $1
}

export PATH=/usr/local/apps/clang/9.0.1/bin:$PATH

if ! [ -x "$(command -v clang-format)" ]; then
  echo 'Error: clang-format is not installed.' >&2
  exit 1
fi

all="yes"
force="no"
dryrun="no"


while test $# -gt 0; do

    # Split --option=value in $opt="--option" and $val="value"

    opt=""
    val=""

    case "$1" in
    --*=*)
      opt=`echo "$1" | sed 's/=.*//'`
      val=`echo "$1" | sed 's/--[_a-zA-Z0-9-]*=//'`
      ;;
    --*)
      opt=$1
      ;;
    *)
      break
      ;;
    esac

    # Parse options
    case "$opt" in
      --help)
        help 0
        exit
        ;;
      --dryrun)
        dryrun="yes"
        ;;
      --force)
        force="yes"
        ;;
      *)
        break
        ;;
    esac
    shift
done

if [ $# -ne 0 ]; then
  all="no"
fi


if ! [[ $(clang-format --version) =~ ${_REQUIRED_CLANG_VERSION} ]]; then
    echo "Error: Require clang-format version: ${_REQUIRED_CLANG_VERSION}"
    echo "    > $(which clang-format) --version"
    echo "      $(clang-format --version)"
    if [[ $force =~ "yes" ]]; then
        echo "--force --> Continue anyway"
    else
        exit 1
    fi
fi

if [[ $all =~ "yes" ]]; then
    echo "Applying $(clang-format --version) to all files ..."

    SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    cd $SCRIPTDIR/../src
    if [[ $dryrun =~ "yes" ]]; then
        echo "+ find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file"
    else
        find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file
    fi
   
    cd $SCRIPTDIR/../atlas_io
    if [[ $dryrun =~ "yes" ]]; then
        echo "+ find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file"
    else
        find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file
    fi

else
    echo "Applying $(clang-format --version) to files ..."
    while test $# -gt 0; do

        if [[ -d $1 ]]; then
            if [[ $dryrun =~ "yes" ]]; then
                echo "+ cd $1"
            fi
            cd $1
            if [[ $dryrun =~ "yes" ]]; then
                echo "+ find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file"
                find . -iname *.h -o -iname *.cc | xargs echo
            else
                find . -iname *.h -o -iname *.cc | xargs clang-format -i -style=file
            fi
            shift
        else
            if [[ $dryrun =~ "yes" ]]; then
                echo "+ clang-format -i -style=file $1"
            else
                clang-format -i -style=file $1
            fi
            shift
        fi
    done
fi

echo "Applying $(clang-format --version) ... done"

