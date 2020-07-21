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


#include "atlas/interpolation/nonlinear/NonLinearFactory.h"

#include <string>

#include "atlas/interpolation/nonlinear/MissingIfAllMissing.h"
#include "atlas/interpolation/nonlinear/MissingIfAnyMissing.h"
#include "atlas/interpolation/nonlinear/MissingIfHeaviestMissing.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


namespace {
void force_link() {
    static struct Link {
        Link() {
            NonLinearFactoryBuilder<MissingIfHeaviestMissing>();
            NonLinearFactoryBuilder<MissingIfAnyMissing>();
            NonLinearFactoryBuilder<MissingIfAllMissing>();
        }
    } link;
}
}  // namespace


const NonLinear* NonLinearFactory::build( const std::string& builder, const util::Config& config ) {
    force_link();
    auto factory = get( builder );
    return factory->make( config );
}


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
