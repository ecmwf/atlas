/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

void throw_OutOfRange(const std::string& class_name, char idx_str, int idx, int max) {
    std::ostringstream msg;
    msg << class_name << " index " << idx << " out of bounds: " << idx << " >= " << max;
    throw_Exception(msg.str(), Here());
}

}  // namespace array
}  // namespace atlas
