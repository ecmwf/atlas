/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


// This file needs to be compiled with an OpenACC capable C compiler that is
// compatible with the Fortran compiler

#ifndef _OPENACC
#error atlas_acc_map_data must be compiled with OpenACC capable compiler
#endif

#include <openacc.h>

void atlas_acc_map_data(void* cpu_ptr, void* gpu_ptr, unsigned long bytes) {
  acc_map_data(cpu_ptr, gpu_ptr, bytes);
}


void atlas_acc_unmap_data(void* cpu_ptr) {
  acc_unmap_data(cpu_ptr);
}

int atlas_acc_is_present(void* cpu_ptr, unsigned long bytes) {
  return acc_is_present(cpu_ptr, bytes);
}

void* atlas_acc_deviceptr(void* cpu_ptr) {
  return acc_deviceptr(cpu_ptr);
}
