// This file needs to be compiled with an OpenACC capable C compiler that is
// compatible with the Fortran compiler

#ifndef _OPENACC
#error atlas_acc_map_data must be compiled with OpenACC capable compiler
#endif

#include <openacc.h>
void atlas_acc_map_data(void* cpu_ptr, void* gpu_ptr, unsigned long size)
{
  acc_map_data(cpu_ptr, gpu_ptr, size);
}
