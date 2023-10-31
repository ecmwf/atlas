/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "RedistributionInterface.h"

//#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

namespace atlas {

// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------

extern "C" {

const Redistribution* atlas__Redistribution__new(
    const functionspace::FunctionSpaceImpl* fspace1, const functionspace::FunctionSpaceImpl* fspace2) {
    return new Redistribution(fspace1, fspace2);
}
}


// ----------------------------------------------------------------------------

}  // namespace atlas
