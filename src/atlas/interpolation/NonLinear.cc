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


#include "atlas/interpolation/NonLinear.h"

#include "atlas/interpolation/nonlinear/NonLinearFactory.h"
#include "atlas/util/Config.h"


namespace atlas {
namespace interpolation {


NonLinear::NonLinear() : Handle( nullptr ) {}


NonLinear::NonLinear( const util::Config& config ) :
    Handle( nonlinear::NonLinearFactory::build( config.getString( "type" ), config ) ) {}


NonLinear::NonLinear( const std::string& type, const util::Config& config ) :
    Handle( nonlinear::NonLinearFactory::build( type, config ) ) {}


}  // namespace interpolation
}  // namespace atlas
