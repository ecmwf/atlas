/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/functionspace/Points.h"

// #include <algorithm>

#include "atlas/array.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/option/Options.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"


namespace atlas {
namespace functionspace {
namespace detail {


Points::Points( const Grid& grid ) :
    lonlat_( "lonlat", array::make_datatype<double>(), array::make_shape( grid.size(), 2 ) ) {
    auto lonlat = array::make_view<double, 2>( lonlat_ );

    idx_t j = 0;
    for ( auto p : grid.lonlat() ) {
        lonlat( j, LON ) = p.lon();
        lonlat( j, LAT ) = p.lat();

        // shouldn't be necessary
        // lonlat( j, LAT ) = std::min(90., std::max(-90., p.lat()));
        ++j;
    }
}


Points::~Points() = default;


Field Points::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | config );
}


Field Points::createField( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) {
        throw_Exception( "datatype missing", Here() );
    }

    Field field( config.getString( "name", "" ), array::DataType( kind ), array::make_shape( size() ) );
    field.set_functionspace( this );
    return field;
}


void Points::haloExchange( const FieldSet&, bool ) const {
    Log::debug() << "Points::haloExchange: ignored" << std::endl;
}


void Points::haloExchange( const Field&, bool ) const {
    Log::debug() << "Points::haloExchange: ignored" << std::endl;
}


idx_t Points::size() const {
    return lonlat_.shape( 0 );
}


size_t Points::footprint() const {
    return sizeof( *this );
}


std::string Points::distribution() const {
    return std::string( "serial" );
}


std::string Points::type() const {
    return "Points";
}


atlas::Field Points::ghost() const {
    if ( not ghost_ ) {
        ghost_ = Field( "ghost", array::make_datatype<int>(), array::make_shape( size() ) );
        array::make_view<int, 1>( ghost_ ).assign( 0 );
    }
    return ghost_;
}


template <typename T>
T Points::IteratorT<T>::operator*() const {
    ATLAS_NOTIMPLEMENTED;
}


template <>
PointXYZ Points::IteratorT<PointXYZ>::operator*() const {
    PointXYZ q;
    Point2 p = {view_( n_, LON ), view_( n_, LAT )};
    util::Earth::convertSphericalToCartesian( p, q );
    return q;
}


template <>
PointLonLat Points::IteratorT<PointLonLat>::operator*() const {
    return {view_( n_, LON ), view_( n_, LAT )};
}


}  // namespace detail


Points::Points( const Grid& grid ) :
    FunctionSpace( new detail::Points( grid ) ),
    functionspace_( dynamic_cast<const detail::Points*>( get() ) ) {}


Points::Points( const FunctionSpace& functionspace ) :
    FunctionSpace( functionspace ),
    functionspace_( dynamic_cast<const detail::Points*>( get() ) ) {}


}  // namespace functionspace
}  // namespace atlas
