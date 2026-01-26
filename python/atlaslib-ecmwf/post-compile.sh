#!/bin/bash

# NOTE we need to strip here because vanilla atlas is for some reason huge
# We don't strip in general since that has caused some linker problem downstream, but here so far seems to work
# NOTE for macos, we additionally fix the install names. For linux we don't need to, since the rpath is already pointing
# to the current directory and there is no strict install name check. The qhull lib was copied already in the precompile
# step

if [ "$(uname)" != "Darwin" ] ; then
    strip --strip-unneeded /tmp/atlas/target/atlas/lib64/libatlas.so
else 
    strip -x /tmp/atlas/target/atlas/lib/libatlas.dylib
    QHULL_NAME=$(otool -l /tmp/atlas/target/atlas/lib/libatlas.dylib | grep qhull | sed 's/.*name \(.*\) (offset.*)/\1/')
    QHULL_BASE=$(basename $QHULL_NAME)
    install_name_tool -change $QHULL_NAME '@rpath/'$QHULL_BASE /tmp/atlas/target/atlas/lib/libatlas.dylib

    # NOTE fftw currently disabled due to licensing reasons
    # FFTW_NAME=$(otool -l /tmp/atlas/target/atlas/lib/libatlas.dylib | grep fftw | sed 's/.*name \(.*\) (offset.*)/\1/')
    # FFTW_BASE=$(basename $FFTW_NAME)
    # install_name_tool -change $FFTW_NAME '@rpath/'$FFTW_BASE /tmp/atlas/target/atlas/lib/libatlas.dylib
fi
