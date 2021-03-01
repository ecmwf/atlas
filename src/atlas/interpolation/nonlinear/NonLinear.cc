/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/nonlinear/NonLinear.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {

void force_link_missing();

const NonLinear* NonLinearFactory::build( const std::string& builder, const NonLinearFactory::Config& config ) {
    force_link_missing();
    return get( builder )->make( config );
}


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
