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


#pragma once

#include "atlas/field/MissingValue.h"
#include "atlas/interpolation/nonlinear/NonLinear.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


struct Missing : NonLinear {
private:
    bool applicable( const Field& f ) const override { return field::MissingValue( f ); }
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
