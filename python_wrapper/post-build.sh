#!/bin/bash

# NOTE auditwheel is problematic since it changes libnames -- all is well from
# the pov # of this very package's libs, but subsequent packages compiled with
# this as a dependency end up not working
# if [ "$(uname)" != "Darwin" ] ; then
#   auditwheel repair -w /tmp/eccodes/auditwheel /tmp/eccodes/build/wheel/*whl
# fi

# TODO macos we actually need to use the `otool` here
