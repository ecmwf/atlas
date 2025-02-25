/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Event.h"

#include <iostream>

#include "pluto/pluto_config.h"
#include "hic/hic.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

void wait(const Event& event) {
    if constexpr (LOG) {
        std::cout << "               = hicEventSynchronize(event:"<<event.value<<")" << std::endl;
    }
}

}
