/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include "pluto/pluto.h"
#include "hic/hic.h"

HIC_HOST_DEVICE
void print() {
  printf("is_on_device = %d\n", int(pluto::is_on_device()));
}

void print_on_host()   { print(); }

HIC_GLOBAL
void print_on_device() { print(); }

int main(int argc, char* argv[]) {
  print_on_host();
  if( pluto::devices() == 0) {
    std::cout << "No devices present" << std::endl;
    return 0;
  }
  #if HIC_COMPILER
    print_on_device<<<1,1>>>();
  #else
    std::cout << "Cannot launch kernel 'print_on_device' as compiler does not support it (HIC_COMPILER=0)." << std::endl;
  #endif
  pluto::wait();
}
