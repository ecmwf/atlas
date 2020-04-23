/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphere.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <numeric>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"
#include "eckit/utils/Translator.h"

#include "atlas/domain/Domain.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"
#include "atlas/util/UnitSphere.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

  /*

     Cubed-sphere Panel Arrangement:

                .......
               |       :
               |   3   :
               |       :
                --->---
      *.......  .......  -------  -------
      |       :|       :|       :|       :
      |   1   :|   2   :v   4   :v   5   :
      |       :|       :|       :|       :
       --->---  --->--- *.......  .......
                                  -------
                                 |       :
                                 v   6   :
                                 |       :
                                  .......

       Key
       Solid lines: left and bottom edges of panel
       Dotted lines: right and top edges of panel
       > < v: direction of increasing index in first dimension
       *: location of two extra points (ngrid = 6 * (NCube+1) * (NCube+1) + 2)

  */

static eckit::Translator<std::string, int> to_int;

static Domain domain( const Grid::Config& grid ) {
    Grid::Config config;
    if ( grid.get( "domain", config ) ) {
        return Domain( config );
    }
    return Domain();
}

std::string CubedSphere::static_type() {
    return "cubedsphere";
}

std::string CubedSphere::name() const {
    return name_;
}

CubedSphere::CubedSphere( int N, Projection p ) :
    CubedSphere( CubedSphere::static_type(), N, p ) {
    }

CubedSphere::CubedSphere( const std::string& name, int CubeNx, Projection projection ) :
    Grid(),
    name_( name ) {
    // Copy members
    CubeNx_ = CubeNx;
    projection_ = projection ? projection : Projection();

    // Total number of nodes
    npts_ = 6 * CubeNx * CubeNx + 2;

    // Domain
    domain_ = computeDomain();

    // x-y starts
    xs_[0] = 0*CubeNx;
    xs_[1] = 1*CubeNx;
    xs_[2] = 1*CubeNx;
    xs_[3] = 2*CubeNx;
    xs_[4] = 3*CubeNx;
    xs_[5] = 3*CubeNx;

    ys_[0] = 1*CubeNx;
    ys_[1] = 1*CubeNx;
    ys_[2] = 2*CubeNx;
    ys_[3] = 2*CubeNx-1;
    ys_[4] = 2*CubeNx-1;
    ys_[5] = 1*CubeNx-1;

    xe_[0] = 1*CubeNx-1;
    xe_[1] = 2*CubeNx-1;
    xe_[2] = 2*CubeNx-1;
    xe_[3] = 3*CubeNx-1;
    xe_[4] = 4*CubeNx-1;
    xe_[5] = 4*CubeNx-1;

    ye_[0] = 2*CubeNx-1;
    ye_[1] = 2*CubeNx-1;
    ye_[2] = 3*CubeNx-1;
    ye_[3] = 1*CubeNx+1;
    ye_[4] = 1*CubeNx+1;
    ye_[5] = 0*CubeNx+1;

    array::ArrayT<int> xArray_( CubeNx+1, CubeNx+1, 6 );
    //array::ArrayView<int, 3> isGhost = array::make_view<int, 3>( isGhostArray );

}

Domain CubedSphere::computeDomain() const {return GlobalDomain(); }

CubedSphere::~CubedSphere() = default;

void CubedSphere::print( std::ostream& os ) const {
    os << "CubedSphere(Name:" << name() << ")";
}

std::string CubedSphere::type() const {
    return static_type();
}

void CubedSphere::hash( eckit::Hash& h ) const {

    h.add("CubedSphere");
    h.add(CubeNx_);

    // also add projection information
    //projection().hash( h );

    // also add domain information, even though already encoded in grid.
    domain().hash( h );

}

RectangularLonLatDomain CubedSphere::lonlatBoundingBox() const {
    return projection_ ? projection_.lonlatBoundingBox( computeDomain() ) : domain();
}

Grid::Spec CubedSphere::spec() const {
    Grid::Spec grid_spec;

    if ( name() == "cubedsphere" ) {
        grid_spec.set( "type", type() );
    }
    else {
        grid_spec.set( "name", name() );
    }
    grid_spec.set( "projection", projection().spec() );
    return grid_spec;
}

// --------------------------------------------------------------------

#if 1
namespace {  // anonymous

//static class cubedsphere : public GridBuilder {
//    using Implementation = atlas::Grid::Implementation;
//    using Config         = Grid::Config;
//
//public:
//    cubedsphere() : GridBuilder( "cubedsphere" ) {}
//
//    void print( std::ostream& os ) const override {
//        os << std::left << std::setw( 20 ) << " "
//           << "CubedSphere grid";
//    }
//
//    const Implementation* create( const std::string& /* name */, const Config& ) const override {
//        throw_NotImplemented( "Cannot create cubedsphere grid from name", Here() );
//    }
//
//    const Implementation* create( const Config& config ) const override {
//        Projection projection;
//
//        Config config_proj;
//        if ( config.get( "projection", config_proj ) ) {
//            projection = Projection( config_proj );
//        }
//
//        int N;
//        config.get( "NumberCubeFaces:", N );
//
//        return new CubedSphereGrid::grid_t( N, projection );
//    }
//
//} cubedsphere_;

}  // anonymous namespace
#endif

// --------------------------------------------------------------------

extern "C" {

idx_t atlas__grid__CubedSphere__ny( CubedSphere* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    return This->ny();
}

idx_t atlas__grid__CubedSphere__nx( CubedSphere* This, idx_t jlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    return This->nx( jlat );
}

void atlas__grid__CubedSphere__nx_array( CubedSphere* This, const idx_t*& nx_array, idx_t& size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    nx_array = This->nx().data();
    size     = idx_t( This->nx().size() );
}

idx_t atlas__grid__CubedSphere__size( CubedSphere* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    return This->size();
}

double atlas__grid__CubedSphere__y( CubedSphere* This, idx_t i, idx_t j, idx_t t ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    return This->y( i, j, t );
}

double atlas__grid__CubedSphere__x( CubedSphere* This, idx_t i, idx_t j, idx_t t ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    return This->x( i, j, t );
}

void atlas__grid__CubedSphere__xy( CubedSphere* This, idx_t i, idx_t j, idx_t t, double crd[] ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    This->xy( i, j, t, crd );
}

void atlas__grid__CubedSphere__lonlat( CubedSphere* This, idx_t i, idx_t j, idx_t t, double crd[] ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    This->lonlat( i, j, t, crd );
}

void atlas__grid__CubedSphere__y_array( CubedSphere* This, const double*& y_array, idx_t& size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    y_array = This->y().data();
    size    = idx_t( This->y().size() );
}

const CubedSphere* atlas__grid__CubedSphere( char* identifier ) {
    const CubedSphere* grid = dynamic_cast<const CubedSphere*>( Grid::create( std::string( identifier ) ) );
    ATLAS_ASSERT( grid != nullptr );
    return grid;
}

const CubedSphere* atlas__grid__CubedSphere__config( util::Config* conf ) {
    ATLAS_ASSERT( conf != nullptr );
    const CubedSphere* grid = dynamic_cast<const CubedSphere*>( Grid::create( *conf ) );
    ATLAS_ASSERT( grid != nullptr );
    return grid;
}

void atlas__grid__CubedSphere__delete( CubedSphere* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_CubedSphereGrid" );
    delete This;
}


}

namespace {
GridFactoryBuilder<CubedSphere> __register_CubedSphere( CubedSphere::static_type() );
}

// -------------------------------------------------------------------------------------------------

// Specialization based on type of projection
// ------------------------------------------

static class cubedsphere_equiangular : public GridBuilder {
public:
    cubedsphere_equiangular() : GridBuilder( "cubedsphere_equiangular", {"^[Cc][Ss][Ee][Aa]([0-9]+)$"}, {"CSEA<cubedsphere>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CSEA<FaceNx>"
           << "Cubed sphere, equiangular";
    }

    const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config gridconf( config );
            int CubeNx = to_int( matches[0] );
            gridconf.set( "type", type() );
            gridconf.set( "CubeNx", CubeNx );
            return create( gridconf );
        }
        return nullptr;
    }

    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int CubeNx;
        config.get( "CubeNx", CubeNx );
        util::Config projconf;
        projconf.set("type", "cubedsphere_equiangular");
        projconf.set("CubeNx", CubeNx);
        return new CubedSphereGrid::grid_t( "CSEA" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_equiangular_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant : public GridBuilder {
public:
    cubedsphere_equidistant() : GridBuilder( "cubedsphere_equidistant", {"^[Cc][Ss][Ee][Dd]([0-9]+)$"}, {"CSED<cubedsphere>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CSED<FaceNx>"
           << "Cubed sphere, equidistant";
    }

    const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config gridconf( config );
            int CubeNx = to_int( matches[0] );
            gridconf.set( "type", type() );
            gridconf.set( "CubeNx", CubeNx );
            return create( gridconf );
        }
        return nullptr;
    }

    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int CubeNx;
        config.get( "CubeNx", CubeNx );
        util::Config projconf;
        projconf.set("type", "cubedsphere_equidistant");
        projconf.set("CubeNx", CubeNx);
        return new CubedSphereGrid::grid_t( "CSED" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_equidistant_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant_fv3 : public GridBuilder {
public:
    cubedsphere_equidistant_fv3() : GridBuilder( "cubedsphere_equidistant_fv3", {"^[Cc][Ss][Ee][Dd][Ff][Vv]3([0-9]+)$"}, {"CSFV3<cubedsphere>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CSEDFV3<FaceNx>"
           << "Cubed sphere, equidistant FV3";
    }

    const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config gridconf( config );
            int CubeNx = to_int( matches[0] );
            gridconf.set( "type", type() );
            gridconf.set( "CubeNx", CubeNx );
            return create( gridconf );
        }
        return nullptr;
    }

    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int CubeNx;
        config.get( "CubeNx", CubeNx );
        util::Config projconf;
        projconf.set("type", "cubedsphere_equidistant_fv3");
        projconf.set("CubeNx", CubeNx);
        return new CubedSphereGrid::grid_t( "CSED" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_equidistant_fv3_;

// -------------------------------------------------------------------------------------------------

void force_link_CubedSphere() {
    cubedsphere_equiangular_.force_link();
    cubedsphere_equidistant_.force_link();
    cubedsphere_equidistant_fv3_.force_link();
  }

// -------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
