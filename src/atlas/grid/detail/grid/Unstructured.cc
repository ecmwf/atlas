/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/grid/Unstructured.h"

#include <limits>
#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "eckit/memory/Builder.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

eckit::ConcreteBuilderT1<Grid, Unstructured> builder_Unstructured( Unstructured::static_type() );

Unstructured::Unstructured( const Mesh& m ) : Grid(), points_( new std::vector<PointXY>( m.nodes().size() ) ) {
    util::Config config_domain;
    config_domain.set( "type", "global" );
    domain_ = Domain( config_domain );

    auto xy                 = array::make_view<double, 2>( m.nodes().xy() );
    std::vector<PointXY>& p = *points_;
    const size_t npts       = p.size();

    for ( size_t n = 0; n < npts; ++n ) {
        p[n].assign( xy( n, XX ), xy( n, YY ) );
    }
}

Unstructured::Unstructured( const util::Config& p ) : Grid() {
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

    for ( size_t i = 0; i < pts.size(); i++ ) {
        const PointXY& p = pts[i];
        h << p.x() << p.y();
    }

    projection().hash( h );
}

size_t Unstructured::size() const {
    ASSERT( points_ );
    return points_->size();
}

Grid::Spec Unstructured::spec() const {
    if ( cached_spec_ ) return *cached_spec_;

    cached_spec_.reset( new Grid::Spec );

    cached_spec_->set( "type", static_type() );

    cached_spec_->set( "domain", domain().spec() );
    cached_spec_->set( "projection", projection().spec() );

    std::unique_ptr<IteratorXY> it( xy_begin() );
    std::vector<double> coords( 2 * size() );
    size_t c( 0 );
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
