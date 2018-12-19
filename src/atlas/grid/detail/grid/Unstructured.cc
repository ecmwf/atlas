/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/grid/Unstructured.h"

#include <initializer_list>
#include <limits>
#include <memory>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/option.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

//static eckit::ConcreteBuilderT1<Grid, Unstructured> builder_Unstructured( Unstructured::static_type() );

namespace {
static GridFactoryBuilder<Unstructured> __register_Unstructured( Unstructured::static_type() );
}


Unstructured::Unstructured( const Mesh& m ) : Grid(), points_( new std::vector<PointXY>( m.nodes().size() ) ) {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );

    auto xy                 = array::make_view<double, 2>( m.nodes().xy() );
    std::vector<PointXY>& p = *points_;
    const idx_t npts        = static_cast<idx_t>( p.size() );

    for ( idx_t n = 0; n < npts; ++n ) {
        p[n].assign( xy( n, XX ), xy( n, YY ) );
    }
}


namespace {
class Normalise {
public:
    Normalise( const RectangularDomain& domain ) :
        degrees_( domain.units() == "degrees" ),
        normalise_( domain.xmin(), domain.xmax() ) {}

    double operator()( double x ) const {
        if ( degrees_ ) { x = normalise_( x ); }
        return x;
    }

private:
    const bool degrees_;
    NormaliseLongitude normalise_;
};
}  // namespace


Unstructured::Unstructured( const Grid& grid, Domain domain ) : Grid() {
    domain_ = domain;
    points_.reset( new std::vector<PointXY> );
    points_->reserve( grid.size() );
    if ( not domain_ ) { domain_ = Domain( option::type( "global" ) ); }
    atlas::grid::IteratorXY it( grid.xy_begin() );
    PointXY p;
    if ( RectangularDomain( domain_ ) ) {
        auto normalise = Normalise( RectangularDomain( domain_ ) );
        while ( it.next( p ) ) {
            p.x() = normalise( p.x() );
            if ( domain_.contains( p ) ) { points_->emplace_back( p ); }
        }
    }
    else if ( ZonalBandDomain( domain_ ) ) {
        while ( it.next( p ) ) {
            if ( domain_.contains( p ) ) { points_->emplace_back( p ); }
        }
    }
    else {
        while ( it.next( p ) ) {
            points_->emplace_back( p );
        }
    }
    points_->shrink_to_fit();
}

Unstructured::Unstructured( const util::Config& ) : Grid() {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );
    NOTIMP;
}

Unstructured::Unstructured( std::vector<PointXY>* pts ) : Grid(), points_( pts ) {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );
}

Unstructured::Unstructured( std::vector<PointXY>&& pts ) :
    Grid(),
    points_( new std::vector<PointXY>( std::move( pts ) ) ) {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );
}

Unstructured::Unstructured( std::initializer_list<PointXY> initializer_list ) :
    Grid(),
    points_( new std::vector<PointXY>( initializer_list ) ) {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );
}

Unstructured::~Unstructured() {}

Grid::uid_t Unstructured::name() const {
    if ( shortName_.empty() ) {
        std::ostringstream s;
        s << "unstructured." << Grid::hash().substr( 0, 7 );
        shortName_ = s.str();
    }
    return shortName_;
}

void Unstructured::hash( eckit::Hash& h ) const {
    ASSERT( points_ );

    const std::vector<PointXY>& pts = *points_;
    h.add( &pts[0], sizeof( PointXY ) * pts.size() );

    for ( idx_t i = 0, N = static_cast<idx_t>( pts.size() ); i < N; i++ ) {
        const PointXY& p = pts[i];
        h << p.x() << p.y();
    }

    projection().hash( h );
}

idx_t Unstructured::size() const {
    ASSERT( points_ );
    return static_cast<idx_t>( points_->size() );
}

Grid::Spec Unstructured::spec() const {
    if ( cached_spec_ ) return *cached_spec_;

    cached_spec_.reset( new Grid::Spec );

    cached_spec_->set( "type", static_type() );

    cached_spec_->set( "domain", domain().spec() );
    cached_spec_->set( "projection", projection().spec() );

    std::unique_ptr<IteratorXY> it( xy_begin() );
    std::vector<double> coords( 2 * size() );
    idx_t c( 0 );
    PointXY xy;
    while ( it->next( xy ) ) {
        coords[c++] = xy.x();
        coords[c++] = xy.y();
    }

    cached_spec_->set( "xy", coords );

    return *cached_spec_;
}

void Unstructured::print( std::ostream& os ) const {
    os << "Unstructured(Npts:" << size() << ")";
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
