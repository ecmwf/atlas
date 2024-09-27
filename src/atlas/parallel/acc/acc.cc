/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "acc.h"

#include "atlas/library/defines.h"

#if ATLAS_HAVE_ACC
#include "pluto/pluto.h"
#include "atlas_acc_support/atlas_acc.h"
#endif

namespace atlas::acc {

int devices() {
#if ATLAS_HAVE_ACC
    static int num_devices = [](){
        if (pluto::devices() == 0) {
            return 0;
        }
        auto devicetype = atlas_acc_get_device_type();
        int _num_devices = atlas_acc_get_num_devices();
        if (_num_devices == 1 && devicetype == atlas_acc_device_host) {
          --_num_devices;
	    }
	    return _num_devices;
    }();
    return num_devices;
#else
    return 0;
#endif
}

void map(void* host_data, void* device_data, std::size_t bytes) {
#if ATLAS_HAVE_ACC
    atlas_acc_map_data(host_data, device_data, bytes);
#endif
}
void unmap(void* host_data) {
#if ATLAS_HAVE_ACC
    atlas_acc_unmap_data(host_data);
#endif
}

bool is_present(void* host_data, std::size_t bytes) {
#if ATLAS_HAVE_ACC
    return atlas_acc_is_present(host_data, bytes);
#else
    return false;
#endif
}

void* deviceptr(void* host_data) {
#if ATLAS_HAVE_ACC
    return atlas_acc_deviceptr(host_data);
#else
    return nullptr;
#endif
}

CompilerId compiler_id() {
#if ATLAS_HAVE_ACC
    static CompilerId id = []() {
        switch (atlas_acc_compiler_id()) {
            case atlas_acc_compiler_id_cray: return CompilerId::cray;
            default: return CompilerId::unknown;
        }
    }();
    return id;
#endif
    return CompilerId::unknown;
}


}

