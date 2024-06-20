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

#ifdef __cplusplus
extern "C" {
#endif

void atlas_acc_map_data(void* cpu_ptr, void* gpu_ptr, unsigned long size);
void atlas_acc_unmap_data(void* cpu_ptr);

#ifdef __cplusplus
}
#endif
