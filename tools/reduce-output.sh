#!/bin/bash

# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


# Abort on Error
set -e

PING_SLEEP=30s

dump_output() {
   echo " ++ Tailing the last 100 lines of output from $BUILD_OUTPUT"
   tail -100 $BUILD_OUTPUT  
}
error_handler() {
  echo ERROR: An error was encountered with the build.
  kill $PING_LOOP_PID
  dump_output
  exit 1
}
# If an error occurs, run our error handler to output a tail of the build
trap 'error_handler' ERR

# Set up a repeating loop to display some output regularly.
bash -c "while true; do sleep $PING_SLEEP; echo \" ++ \$(date) - running ... \"; done" &
PING_LOOP_PID=$!
BUILD_OUTPUT=build-$PING_LOOP_PID.out
touch $BUILD_OUTPUT
echo " + $@"
echo " ++ Output redirected to $BUILD_OUTPUT"
$@ >> $BUILD_OUTPUT 2>&1

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
