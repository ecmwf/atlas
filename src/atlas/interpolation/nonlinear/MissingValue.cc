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


#include "atlas/interpolation/nonlinear/MissingValue.h"

#include <cmath>

#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


namespace {
double config_value( const MissingValue::Config& c ) {
    double value;
    ATLAS_ASSERT( c.get( "missing_value", value ) );
    return value;
}


double config_epsilon( const MissingValue::Config& c ) {
    double value = 0.;
    c.get( "missing_value_epsilon", value );
    return value;
}
}  // namespace


/// @brief Missing value if NaN
struct MissingValueNaN : MissingValue {
    MissingValueNaN( const Config& ) {}
    bool operator()( const double& value ) const override { return std::isnan( value ); }
    bool isnan() const override { return true; }
    static std::string static_type() { return "nan"; }
};


/// @brief Missing value if comparing equally to pre-defined value
struct MissingValueEquals : MissingValue {
    MissingValueEquals( const Config& config ) :
        missingValue_( config_value( config ) ), missingValue2_( missingValue_ ) {
        ATLAS_ASSERT( missingValue_ == missingValue2_ );  // this succeeds
    }

    MissingValueEquals( double missingValue ) : missingValue_( missingValue ), missingValue2_( missingValue_ ) {
        ATLAS_ASSERT( !std::isnan( missingValue2_ ) );
    }

    bool operator()( const double& value ) const override {
        // ATLAS_ASSERT(missingValue_ == missingValue2_);  // this fails when not using missingValue2_ (copy ellision
        // problem on POD!?)
        return value == missingValue2_;
    }

    bool isnan() const override { return false; }

    static std::string static_type() { return "equals"; }

    const double missingValue_;
    const double missingValue2_;
};


/// @brief Missing value if comparing approximately to pre-defined value
struct MissingValueApprox : MissingValue {
    MissingValueApprox( const Config& config ) :
        missingValue_( config_value( config ) ), epsilon_( config_epsilon( config ) ) {}
    MissingValueApprox( double missingValue, double epsilon ) : missingValue_( missingValue ), epsilon_( epsilon ) {
        ATLAS_ASSERT( !std::isnan( missingValue_ ) );
        ATLAS_ASSERT( epsilon_ >= 0. );
    }

    bool operator()( const double& value ) const override {
        return eckit::types::is_approximately_equal( value, missingValue_, epsilon_ );
    }

    bool isnan() const override { return false; }

    static std::string static_type() { return "approximately-equals"; }

    const double missingValue_;
    const double epsilon_;
};


MissingValue::~MissingValue() = default;


namespace {
MissingValueFactoryBuilder<MissingValueNaN> __mv1;
MissingValueFactoryBuilder<MissingValueApprox> __mv2;
MissingValueFactoryBuilder<MissingValueEquals> __mv3;
}  // namespace


const MissingValue* MissingValueFactory::build( const std::string& builder, const Config& config ) {
    return has( builder ) ? get( builder )->make( config ) : nullptr;
}


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
