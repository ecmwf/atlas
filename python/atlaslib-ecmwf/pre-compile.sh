#!/bin/bash

# the procedure for adding a new ext dependency to be bundled in here:
# - add git checkout, compile, etc
# - ensure the version ends up in python/atlaslib-ecmwf/src/versions.txt
# - ensure the licence ends up in python/atlaslib-ecmwf/src/copying/, and fname is referenced in copying/list.json
# - ensure the .so ends up in target/lib64/ with the expected libname
# - validate that the resulting wheel contains all the above
# additionally, make sure this script is aligned with /buildscripts/compile.sh and /buildscripts/wheel-linux.sh,
# in particular when it comes to install targets and package data, etc

set -euo pipefail

mkdir -p python/atlaslib-ecmwf/src/copying

# NOTE we dont use that since we dont want to work with the brew installation as it is unreliable
# bash ./tools/install-qhull.sh --prefix /tmp/qhull/target
mkdir -p /tmp/qhull/build && cd /tmp/qhull/build
git clone https://github.com/qhull/qhull /tmp/qhull/src
cmake /tmp/qhull/src -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/tmp/qhull/target
make -j8 && make install
cd -

# copy the libs, instead of having auditwheel done it later. This is a bit risky because cmake will later write in this
# very directory... but it works
if [ "$(uname)" != "Darwin" ] ; then
    mkdir -p /tmp/atlas/target/atlas/lib64
    cp /tmp/qhull/target/lib/libqhull_r.* /tmp/atlas/target/atlas/lib64/
else
    mkdir -p /tmp/atlas/target/atlas/lib
    cp /tmp/qhull/target/lib/libqhull_r.* /tmp/atlas/target/atlas/lib/
fi

wget https://raw.githubusercontent.com/qhull/qhull/master/COPYING.txt -O python/atlaslib-ecmwf/src/copying/libqhull.txt
echo '{"libqhull_r": {"path": "copying/libqhull.txt", "home": "https://github.com/qhull/qhull"}}' > python/atlaslib-ecmwf/src/copying/list.json
