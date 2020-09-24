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


#include "atlas/field/detail/MissingValue.h"

#include <cmath>
#include <type_traits>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/DataType.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"


namespace atlas {
namespace field {
namespace detail {


namespace {
static const std::string type_key    = "missing_value_type";
static const std::string value_key   = "missing_value";
static const std::string epsilon_key = "missing_value_epsilon";


template <typename T>
T config_value( const MissingValue::Config& c ) {
    T value;
    ATLAS_ASSERT( c.get( value_key, value ) );
    return value;
}


template <typename T>
T config_epsilon( const MissingValue::Config& c ) {
    T value = 0.;
    c.get( epsilon_key, value );
    return value;
}
}  // namespace


/**
 * @brief Missing value if NaN
 */
template <typename T>
struct MissingValueNaN : MissingValue {
    MissingValueNaN( const Config& ) { ATLAS_ASSERT( std::is_floating_point<T>::value ); }
    bool operator()( const T& value ) const override { return std::isnan( value ); }
    bool isnan() const override { return true; }
    void metadata( Field& field ) const override { field.metadata().set( type_key, static_type() ); }
    static std::string static_type() { return "nan"; }
};


/**
 * @brief Missing value if comparing equally to pre-defined value
 */
template <typename T>
struct MissingValueEquals : MissingValue {
    MissingValueEquals( const Config& config ) : MissingValueEquals( config_value<T>( config ) ) {}

    MissingValueEquals( T missingValue ) : missingValue_( missingValue ), missingValue2_( missingValue_ ) {
        ATLAS_ASSERT( missingValue_ == missingValue2_ );  // FIXME this succeeds
        ATLAS_ASSERT( !std::isnan( missingValue2_ ) );
    }

    bool operator()( const T& value ) const override {
        // ATLAS_ASSERT(missingValue_ == missingValue2_);  // FIXME this fails (copy ellision problem on POD!?)
        return value == missingValue2_;
    }

    bool isnan() const override { return false; }

    void metadata( Field& field ) const override {
        field.metadata().set( type_key, static_type() );
        field.metadata().set( value_key, missingValue2_ );
    }

    static std::string static_type() { return "equals"; }

    const T missingValue_;
    const T missingValue2_;
};


/**
 * @brief Missing value if comparing approximately to pre-defined value
 */
template <typename T>
struct MissingValueApprox : MissingValue {
    MissingValueApprox( const Config& config ) :
        MissingValueApprox( config_value<T>( config ), config_epsilon<T>( config ) ) {}

    MissingValueApprox( T missingValue, T epsilon ) : missingValue_( missingValue ), epsilon_( epsilon ) {
        ATLAS_ASSERT( !std::isnan( missingValue_ ) );
        ATLAS_ASSERT( std::is_floating_point<T>::value );
        ATLAS_ASSERT( epsilon_ >= 0. );
    }

    bool operator()( const T& value ) const override {
        return eckit::types::is_approximately_equal( value, missingValue_, epsilon_ );
    }

    bool isnan() const override { return false; }

    void metadata( Field& field ) const override {
        field.metadata().set( type_key, static_type() );
        field.metadata().set( value_key, missingValue_ );
        field.metadata().set( epsilon_key, epsilon_ );
    }

    static std::string static_type() { return "approximately-equals"; }

    const T missingValue_;
    const T epsilon_;
};


const MissingValue* MissingValueFactory::build( const std::string& builder, const Config& config ) {
    return has( builder ) ? get( builder )->make( config ) : nullptr;
}


#define B MissingValueFactoryBuilder
#define M1 MissingValueNaN
#define M2 MissingValueEquals
#define M3 MissingValueApprox
#define T array::DataType::str


static B<M1<double>> __mv1( M1<double>::static_type() );
static B<M1<double>> __mv2( M1<double>::static_type() + "-" + T<double>() );
static B<M1<float>> __mv3( M1<float>::static_type() + "-" + T<float>() );

static B<M2<double>> __mv4( M2<double>::static_type() );
static B<M2<double>> __mv5( M2<double>::static_type() + "-" + T<double>() );
static B<M2<float>> __mv6( M2<float>::static_type() + "-" + T<float>() );
static B<M2<int>> __mv7( M2<int>::static_type() + "-" + T<int>() );
static B<M2<long>> __mv8( M2<long>::static_type() + "-" + T<long>() );
static B<M2<unsigned long>> __mv9( M2<unsigned long>::static_type() + "-" + T<unsigned long>() );

static B<M3<double>> __mv10( M3<double>::static_type() );
static B<M3<double>> __mv11( M3<double>::static_type() + "-" + T<double>() );
static B<M3<float>> __mv12( M3<float>::static_type() + "-" + T<float>() );


#undef T
#undef M3
#undef M2
#undef M1
#undef B


}  // namespace detail
}  // namespace field
}  // namespace atlas
