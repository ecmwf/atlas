#!/bin/bash
if [ "$(uname)" != "Darwin" ] ; then
    strip --strip-unneeded /tmp/atlas/target/atlas/lib64/libatlas.so
else 
    strip --strip-unneeded /tmp/atlas/target/atlas/lib/libatlas.so # TODO test this
fi
