/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    atlas_acc_device_host = 0,
    atlas_acc_device_not_host = 1
} atlas_acc_device_t;

void atlas_acc_map_data(void* cpu_ptr, void* gpu_ptr, size_t bytes);
void atlas_acc_unmap_data(void* cpu_ptr);
int atlas_acc_is_present(void* cpu_ptr, size_t bytes);
void* atlas_acc_deviceptr(void* cpu_ptr);
atlas_acc_device_t atlas_acc_get_device_type();
int atlas_acc_get_num_devices();

typedef enum {
    atlas_acc_compiler_id_unknown = 0,
    atlas_acc_compiler_id_nvidia  = 1,
    atlas_acc_compiler_id_cray    = 2
} atlas_acc_compiler_id_t;

atlas_acc_compiler_id_t atlas_acc_compiler_id();

#ifdef __cplusplus
}
#endif
