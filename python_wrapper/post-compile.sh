#!/bin/bash

# NOTE we need to strip here because vanilla atlas is for some reason huge
# We don't strip in general since that has caused some linker problem downstream, but here so far seems to work

if [ "$(uname)" != "Darwin" ] ; then
    strip --strip-unneeded /tmp/atlas/target/atlas/lib64/libatlas.so
else 
    strip --strip-unneeded -x /tmp/atlas/target/atlas/lib/libatlas.dylib
fi
