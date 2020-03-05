/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/functionspace/PointCloud.h"
#include "atlas/array.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/option/Options.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace functionspace {

namespace detail {

PointCloud::PointCloud( PointXY, const std::vector<PointXY>& points ) : PointCloud( points ) {}

PointCloud::PointCloud( const std::vector<PointXY>& points ) {
    lonlat_     = Field( "lonlat", array::make_datatype<double>(), array::make_shape( points.size(), 2 ) );
    auto lonlat = array::make_view<double, 2>( lonlat_ );
    for ( idx_t j = 0, size = points.size(); j < size; ++j ) {
        lonlat( j, 0 ) = points[j].x();
        lonlat( j, 1 ) = points[j].y();
    }
}

PointCloud::PointCloud( PointXYZ, const std::vector<PointXYZ>& points ) {
    lonlat_       = Field( "lonlat", array::make_datatype<double>(), array::make_shape( points.size(), 2 ) );
    vertical_     = Field( "vertical", array::make_datatype<double>(), array::make_shape( points.size() ) );
    auto lonlat   = array::make_view<double, 2>( lonlat_ );
    auto vertical = array::make_view<double, 1>( vertical_ );
    for ( idx_t j = 0, size = points.size(); j < size; ++j ) {
        lonlat( j, 0 ) = points[j].x();
        lonlat( j, 1 ) = points[j].y();
        vertical( j )  = points[j].z();
    }
}

PointCloud::PointCloud( const Field& lonlat ) : lonlat_( lonlat ) {}

PointCloud::PointCloud( const Field& lonlat, const Field& ghost ) : lonlat_( lonlat ), ghost_( ghost ) {}

PointCloud::PointCloud( const Grid& grid ) {
    lonlat_     = Field( "lonlat", array::make_datatype<double>(), array::make_shape( grid.size(), 2 ) );
    auto lonlat = array::make_view<double, 2>( lonlat_ );

    idx_t j{0};
    for ( auto p : grid.lonlat() ) {
        lonlat( j, 0 ) = p.lon();
        lonlat( j, 1 ) = p.lat();
        ++j;
    }
}

const Field& PointCloud::ghost() const {
    if ( not ghost_ ) {
        ghost_ = Field( "ghost", array::make_datatype<int>(), array::make_shape( size() ) );
        array::make_view<int, 1>( ghost_ ).assign( 0 );
    }
    return ghost_;
}

Field PointCloud::createField( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) {
        throw_Exception( "datatype missing", Here() );
    }
    auto datatype = array::DataType( kind );

    std::string name;
    config.get( "name", name );
    idx_t levels = levels_;
    config.get( "levels", levels );
    Field field;
    if ( levels ) {
        field = Field( name, datatype, array::make_shape( size(), levels ) );
        field.set_levels( levels );
    }
    else {
        field = Field( name, datatype, array::make_shape( size() ) );
    }
    return field;
}

Field PointCloud::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | config );
}

std::string PointCloud::distribution() const {
    return std::string( "serial" );
}

atlas::functionspace::detail::PointCloud::IteratorXYZ::IteratorXYZ( const atlas::functionspace::detail::PointCloud& fs,
                                                                    bool begin ) :
    fs_( fs ),
    xy_( array::make_view<double, 2>( fs_.lonlat() ) ),
    z_( array::make_view<double, 1>( fs_.vertical() ) ),
    n_( begin ? 0 : fs_.size() ) {}

bool atlas::functionspace::detail::PointCloud::IteratorXYZ::next( PointXYZ& xyz ) {
    if ( n_ < fs_.size() ) {
        xyz.x() = xy_( n_, 0 );
        xyz.y() = xy_( n_, 1 );
        xyz.z() = z_( n_ );
        ++n_;
        return true;
    }
    return false;
}

atlas::functionspace::detail::PointCloud::IteratorXY::IteratorXY( const atlas::functionspace::detail::PointCloud& fs,
                                                                  bool begin ) :
    fs_( fs ), xy_( array::make_view<double, 2>( fs_.lonlat() ) ), n_( begin ? 0 : fs_.size() ) {}

bool atlas::functionspace::detail::PointCloud::IteratorXY::next( PointXY& xyz ) {
    if ( n_ < fs_.size() ) {
        xyz.x() = xy_( n_, 0 );
        xyz.y() = xy_( n_, 1 );
        ++n_;
        return true;
    }
    return false;
}

const PointXY atlas::functionspace::detail::PointCloud::IteratorXY::operator*() const {
    PointXY xy;
    xy.x() = xy_( n_, 0 );
    xy.y() = xy_( n_, 1 );
    return xy;
}

const PointXYZ atlas::functionspace::detail::PointCloud::IteratorXYZ::operator*() const {
    PointXYZ xyz;
    xyz.x() = xy_( n_, 0 );
    xyz.y() = xy_( n_, 1 );
    xyz.z() = z_( n_ );
    return xyz;
}

}  // namespace detail

PointCloud::PointCloud( const FunctionSpace& functionspace ) :
    FunctionSpace( functionspace ), functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}

PointCloud::PointCloud( const Field& points ) :
    FunctionSpace( new detail::PointCloud( points ) ),
    functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}

PointCloud::PointCloud( const std::vector<PointXY>& points ) :
    FunctionSpace( new detail::PointCloud( points ) ),
    functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}

PointCloud::PointCloud( PointXY p, const std::vector<PointXY>& points ) :
    FunctionSpace( new detail::PointCloud( p, points ) ),
    functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}

PointCloud::PointCloud( PointXYZ p, const std::vector<PointXYZ>& points ) :
    FunctionSpace( new detail::PointCloud( p, points ) ),
    functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}

PointCloud::PointCloud( const Grid& grid ) :
    FunctionSpace( new detail::PointCloud( grid ) ),
    functionspace_( dynamic_cast<const detail::PointCloud*>( get() ) ) {}


}  // namespace functionspace
}  // namespace atlas
