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

#include "atlas/interpolation/nonlinear/MissingIfAllMissing.h"
#include "atlas/interpolation/nonlinear/MissingIfAnyMissing.h"
#include "atlas/interpolation/nonlinear/MissingIfHeaviestMissing.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


#define B NonLinearFactoryBuilder
#define M1 MissingIfAllMissing
#define M2 MissingIfAnyMissing
#define M3 MissingIfHeaviestMissing
#define T array::DataType::str


static B<M1<double>> __nl1( M1<double>::static_type() );
static B<M1<double>> __nl2( M1<double>::static_type() + "-" + T<double>() );
static B<M1<float>> __nl3( M1<float>::static_type() + "-" + T<float>() );
static B<M1<int>> __nl4( M1<int>::static_type() + "-" + T<int>() );
static B<M1<long>> __nl5( M1<long>::static_type() + "-" + T<long>() );
static B<M1<unsigned long>> __nl6( M1<unsigned long>::static_type() + "-" + T<unsigned long>() );

static B<M2<double>> __nl7( M2<double>::static_type() );
static B<M2<double>> __nl8( M2<double>::static_type() + "-" + T<double>() );
static B<M2<float>> __nl9( M2<float>::static_type() + "-" + T<float>() );
static B<M2<int>> __nl10( M2<int>::static_type() + "-" + T<int>() );
static B<M2<long>> __nl11( M2<long>::static_type() + "-" + T<long>() );
static B<M2<unsigned long>> __nl12( M2<unsigned long>::static_type() + "-" + T<unsigned long>() );

static B<M3<double>> __nl13( M3<double>::static_type() );
static B<M3<double>> __nl14( M3<double>::static_type() + "-" + T<double>() );
static B<M3<float>> __nl15( M3<float>::static_type() + "-" + T<float>() );
static B<M3<int>> __nl16( M3<int>::static_type() + "-" + T<int>() );
static B<M3<long>> __nl17( M3<long>::static_type() + "-" + T<long>() );
static B<M3<unsigned long>> __nl18( M3<unsigned long>::static_type() + "-" + T<unsigned long>() );


#undef T
#undef M3
#undef M2
#undef M1
#undef B


const NonLinear* NonLinearFactory::build( const std::string& builder, const Config& config ) {
    return get( builder )->make( config );
}


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
