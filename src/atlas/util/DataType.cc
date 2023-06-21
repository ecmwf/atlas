/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DataType.h"

#include <sstream>

#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

void DataType::throw_not_recognised(kind_t kind) {
    std::stringstream msg;
    msg << "kind " << kind << " not recognised.";
    throw_Exception(msg.str(), Here());
}

void DataType::throw_not_recognised(std::string datatype) {
    std::stringstream msg;
    msg << "datatype " << datatype << " not recognised.";
    throw_Exception(msg.str(), Here());
}


//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
