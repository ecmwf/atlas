#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void atlas_acc_map_data( void* cpu_ptr, void* gpu_ptr, unsigned long size );

#ifdef __cplusplus
}
#endif
