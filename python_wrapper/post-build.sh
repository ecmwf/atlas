#!/bin/bash

# NOTE auditwheel is problematic since it changes libnames -- all is well from
# the pov of this very package's libs, but subsequent packages compiled with
# this as a dependency end up not working. Additionally, auditwheel would have
# included eckitlib, which we don't want. Instead we solve the problem in pre-compile
# if [ "$(uname)" != "Darwin" ] ; then
#   auditwheel repair -w /tmp/eccodes/auditwheel /tmp/eccodes/build/wheel/*whl
# fi

# NOTE on macos, we have fixed the problem already in the post-build.sh step
