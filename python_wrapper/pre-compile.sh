#!/bin/bash

# the procedure for adding a new ext dependency to be bundled in here:
# - add git checkout, compile, etc
# - ensure the version ends up in python_wrapper/src/versions.txt
# - ensure the licence ends up in python_wrapper/src/copying/, and fname is referenced in copying/list.json
# - ensure the .so ends up in target/lib64/ with the expected libname
# - validate that the resulting wheel contains all the above
# additionally, make sure this script is aligned with /buildscripts/compile.sh and /buildscripts/wheel-linux.sh,
# in particular when it comes to install targets and package data, etc

set -euo pipefail

mkdir -p python_wrapper/src/copying

bash ./tools/install-qhull.sh --prefix /tmp/qhull/target
# copy the libs, instead of having auditwheel done it later. This is a bit risky because cmake will later write in this
# very directory... but it works
if [ "$(uname)" != "Darwin" ] ; then
    mkdir -p /tmp/atlas/target/atlas/lib64
    cp /tmp/qhull/target/lib/libqhull_r.* /tmp/atlas/target/atlas/lib64/
else
    mkdir -p /tmp/atlas/target/atlas/lib
    for e in $(brew list qhull | grep '/lib/' | grep dylib) ; do
        cp $e /tmp/atlas/target/atlas/lib
    done
fi

wget https://raw.githubusercontent.com/qhull/qhull/master/COPYING.txt -O python_wrapper/src/copying/libqhull.txt
echo '{"libqhull_r": {"path": "copying/libqhull.txt", "home": "https://github.com/qhull/qhull"}}' > python_wrapper/src/copying/list.json



