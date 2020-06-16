/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/StructuredColumns.h"

#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>

#include "eckit/utils/MD5.h"

#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/domain.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/StructuredPartitionPolygon.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/fill.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/detail/Cache.h"

#define REMOTE_IDX_BASE 0

namespace atlas {
namespace functionspace {
namespace detail {

namespace {


template <typename T, typename Field>
array::LocalView<T, 3> make_leveled_view( Field& field ) {
    using namespace array;
    if ( field.levels() ) {
        if ( field.variables() ) {
            return make_view<T, 3>( field ).slice( Range::all(), Range::all(), Range::all() );
        }
        else {
            return make_view<T, 2>( field ).slice( Range::all(), Range::all(), Range::dummy() );
        }
    }
    else {
        if ( field.variables() ) {
            return make_view<T, 2>( field ).slice( Range::all(), Range::dummy(), Range::all() );
        }
        else {
            return make_view<T, 1>( field ).slice( Range::all(), Range::dummy(), Range::dummy() );
        }
    }
}

template <typename T>
std::string checksum_3d_field( const parallel::Checksum& checksum, const Field& field ) {
    auto values = make_leveled_view<const T>( field );
    array::ArrayT<T> surface_field( values.shape( 0 ), values.shape( 2 ) );
    auto surface     = array::make_view<T, 2>( surface_field );
    const idx_t npts = values.shape( 0 );
    atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t j = 0; j < surface.shape( 1 ); ++j ) {
            surface( n, j ) = 0.;
            for ( idx_t l = 0; l < values.shape( 1 ); ++l ) {
                surface( n, j ) += values( n, l, j );
            }
        }
    }
    return checksum.execute( surface.data(), surface_field.stride( 0 ) );
}

}  // namespace


class StructuredColumnsHaloExchangeCache : public util::Cache<std::string, parallel::HaloExchange>,
                                           public grid::detail::grid::GridObserver {
private:
    using Base = util::Cache<std::string, parallel::HaloExchange>;
    StructuredColumnsHaloExchangeCache() : Base( "StructuredColumnsHaloExchangeCache" ) {}

public:
    static StructuredColumnsHaloExchangeCache& instance() {
        static StructuredColumnsHaloExchangeCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& funcspace ) {
        registerGrid( *funcspace.grid().get() );

        creator_type creator = std::bind( &StructuredColumnsHaloExchangeCache::create, &funcspace );
        return Base::get_or_create( key( funcspace ), remove_key( funcspace ), creator );
    }
    void onGridDestruction( grid::detail::grid::Grid& grid ) override { remove( remove_key( grid ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& funcspace ) {
        std::ostringstream key;
        key << "grid[address=" << funcspace.grid().get() << ",halo=" << funcspace.halo()
            << ",periodic_points=" << std::boolalpha << funcspace.periodic_points_
            << ",distribution=" << funcspace.distribution() << "]";
        return key.str();
    }

    static Base::key_type remove_key( const detail::StructuredColumns& funcspace ) {
        return remove_key( *funcspace.grid().get() );
    }

    static Base::key_type remove_key( const grid::detail::grid::Grid& grid ) {
        std::ostringstream key;
        key << "grid[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* funcspace ) {
        value_type* value = new value_type();

        value->setup( array::make_view<int, 1>( funcspace->partition() ).data(),
                      array::make_view<idx_t, 1>( funcspace->remote_index() ).data(), REMOTE_IDX_BASE,
                      funcspace->sizeHalo(), funcspace->sizeOwned() );
        return value;
    }
    ~StructuredColumnsHaloExchangeCache() override = default;
};

class StructuredColumnsGatherScatterCache : public util::Cache<std::string, parallel::GatherScatter>,
                                            public grid::detail::grid::GridObserver {
private:
    using Base = util::Cache<std::string, parallel::GatherScatter>;
    StructuredColumnsGatherScatterCache() : Base( "StructuredColumnsGatherScatterCache" ) {}

public:
    static StructuredColumnsGatherScatterCache& instance() {
        static StructuredColumnsGatherScatterCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& funcspace ) {
        registerGrid( *funcspace.grid().get() );
        creator_type creator = std::bind( &StructuredColumnsGatherScatterCache::create, &funcspace );
        return Base::get_or_create( key( funcspace ), remove_key( funcspace ), creator );
    }
    void onGridDestruction( grid::detail::grid::Grid& grid ) override { remove( remove_key( grid ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& funcspace ) {
        std::ostringstream key;
        key << "grid[address=" << funcspace.grid().get() << ",halo=" << funcspace.halo()
            << ",periodic_points=" << std::boolalpha << funcspace.periodic_points_
            << ",distribution=" << funcspace.distribution() << "]";
        return key.str();
    }

    static Base::key_type remove_key( const detail::StructuredColumns& funcspace ) {
        return remove_key( *funcspace.grid().get() );
    }

    static Base::key_type remove_key( const grid::detail::grid::Grid& grid ) {
        std::ostringstream key;
        key << "grid[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* funcspace ) {
        value_type* value = new value_type();

        value->setup( array::make_view<int, 1>( funcspace->partition() ).data(),
                      array::make_view<idx_t, 1>( funcspace->remote_index() ).data(), REMOTE_IDX_BASE,
                      array::make_view<gidx_t, 1>( funcspace->global_index() ).data(), funcspace->sizeOwned() );
        return value;
    }
    ~StructuredColumnsGatherScatterCache() override = default;
};

class StructuredColumnsChecksumCache : public util::Cache<std::string, parallel::Checksum>,
                                       public grid::detail::grid::GridObserver {
private:
    using Base = util::Cache<std::string, parallel::Checksum>;
    StructuredColumnsChecksumCache() : Base( "StructuredColumnsChecksumCache" ) {}

public:
    static StructuredColumnsChecksumCache& instance() {
        static StructuredColumnsChecksumCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& funcspace ) {
        registerGrid( *funcspace.grid().get() );
        creator_type creator = std::bind( &StructuredColumnsChecksumCache::create, &funcspace );
        return Base::get_or_create( key( funcspace ), remove_key( funcspace ), creator );
    }
    void onGridDestruction( grid::detail::grid::Grid& grid ) override { remove( remove_key( grid ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& funcspace ) {
        std::ostringstream key;
        key << "grid[address=" << funcspace.grid().get() << ",halo=" << funcspace.halo()
            << ",periodic_points=" << std::boolalpha << funcspace.periodic_points_
            << ",distribution=" << funcspace.distribution() << "]";
        return key.str();
    }

    static Base::key_type remove_key( const detail::StructuredColumns& funcspace ) {
        return remove_key( *funcspace.grid().get() );
    }

    static Base::key_type remove_key( const grid::detail::grid::Grid& grid ) {
        std::ostringstream key;
        key << "grid[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* funcspace ) {
        value_type* value = new value_type();
        util::ObjectHandle<parallel::GatherScatter> gather(
            StructuredColumnsGatherScatterCache::instance().get_or_create( *funcspace ) );
        value->setup( gather );
        return value;
    }
    ~StructuredColumnsChecksumCache() override = default;
};


const parallel::GatherScatter& StructuredColumns::gather() const {
    if ( gather_scatter_ ) {
        return *gather_scatter_;
    }
    gather_scatter_ = StructuredColumnsGatherScatterCache::instance().get_or_create( *this );
    return *gather_scatter_;
}

const parallel::GatherScatter& StructuredColumns::scatter() const {
    if ( gather_scatter_ ) {
        return *gather_scatter_;
    }
    gather_scatter_ = StructuredColumnsGatherScatterCache::instance().get_or_create( *this );
    return *gather_scatter_;
}

const parallel::Checksum& StructuredColumns::checksum() const {
    if ( checksum_ ) {
        return *checksum_;
    }
    checksum_ = StructuredColumnsChecksumCache::instance().get_or_create( *this );
    return *checksum_;
}

const parallel::HaloExchange& StructuredColumns::halo_exchange() const {
    if ( halo_exchange_ ) {
        return *halo_exchange_;
    }
    halo_exchange_ = StructuredColumnsHaloExchangeCache::instance().get_or_create( *this );
    return *halo_exchange_;
}

void StructuredColumns::set_field_metadata( const eckit::Configuration& config, Field& field ) const {
    field.set_functionspace( this );

    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            field.metadata().set( "owner", owner );
        }
    }
    field.metadata().set( "global", global );

    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    field.set_levels( levels );

    idx_t variables( 0 );
    config.get( "variables", variables );
    field.set_variables( variables );

    if ( config.has( "type" ) ) {
        field.metadata().set( "type", config.getString( "type" ) );
    }
}

array::DataType StructuredColumns::config_datatype( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) {
        throw_Exception( "datatype missing", Here() );
    }
    return array::DataType( kind );
}

std::string StructuredColumns::config_name( const eckit::Configuration& config ) const {
    std::string name;
    config.get( "name", name );
    return name;
}

idx_t StructuredColumns::config_levels( const eckit::Configuration& config ) const {
    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    return levels;
}

array::ArrayShape StructuredColumns::config_shape( const eckit::Configuration& config ) const {
    array::ArrayShape shape;

    shape.emplace_back( config_size( config ) );

    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    if ( levels > 0 ) {
        shape.emplace_back( levels );
    }

    idx_t variables( 0 );
    config.get( "variables", variables );
    if ( variables > 0 ) {
        shape.emplace_back( variables );
    }

    return shape;
}
size_t StructuredColumns::Map2to1::footprint() const {
    size_t size = sizeof( *this );
    size += data_.size() * sizeof( decltype( data_ )::value_type );
    return size;
}

void StructuredColumns::Map2to1::print( std::ostream& out ) const {
    for ( idx_t j = j_min_; j <= j_max_; ++j ) {
        out << std::setw( 4 ) << j << " : ";
        for ( idx_t i = i_min_; i <= i_max_; ++i ) {
            idx_t v = operator()( i, j );
            if ( v == missing() ) {
                out << std::setw( 6 ) << "X";
            }
            else {
                out << std::setw( 6 ) << v;
            }
        }
        out << '\n';
    }
}

void StructuredColumns::IndexRange::print( std::ostream& out ) const {
    for ( idx_t i = min_; i <= max_; ++i ) {
        idx_t v = operator()( i );
        if ( v == missing() ) {
            out << std::setw( 4 ) << "X";
        }
        else {
            out << std::setw( 4 ) << v;
        }
    }
    out << '\n';
}

idx_t StructuredColumns::config_size( const eckit::Configuration& config ) const {
    idx_t size = size_halo_;
    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            size = ( mpi::rank() == owner ? grid_->size() : 0 );
        }
    }
    return size;
}

std::string StructuredColumns::distribution() const {
    return distribution_;
}

void StructuredColumns::throw_outofbounds( idx_t i, idx_t j ) const {
    std::stringstream ss;
    if ( j < j_begin_halo() || j >= j_end_halo() ) {
        ss << "OutofBounds: j out of range! : (i,j) = (" << i << "," << j << ") --- Expected: " << j_begin_halo()
           << " <= j < " << j_end_halo();
    }
    if ( i < i_begin_halo( j ) || i >= i_end_halo( j ) ) {
        ss << "OutofBounds: i out of range! : (i,j) = (" << i << "," << j << ") --- Expected: " << i_begin_halo( j )
           << " <= i < " << i_end_halo( j );
    }
    throw_Exception( ss.str(), Here() );
}

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
StructuredColumns::StructuredColumns( const Grid& grid, const eckit::Configuration& config ) :
    StructuredColumns::StructuredColumns( grid, grid::Partitioner(), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Partitioner& p,
                                      const eckit::Configuration& config ) :
    StructuredColumns( grid, Vertical( config ), p, config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Distribution& distribution,
                                      const eckit::Configuration& config ) :
    StructuredColumns( grid, distribution, Vertical( config ), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Distribution& distribution,
                                      const Vertical& vertical, const eckit::Configuration& config ) :
    vertical_( vertical ), nb_levels_( vertical_.size() ), grid_( new StructuredGrid( grid ) ) {
    setup( distribution, config );
}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const eckit::Configuration& config ) :
    StructuredColumns( grid, vertical, grid::Partitioner(), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const grid::Partitioner& p,
                                      const eckit::Configuration& config ) :
    vertical_( vertical ), nb_levels_( vertical_.size() ), grid_( new StructuredGrid( grid ) ) {
    ATLAS_TRACE( "StructuredColumns constructor" );

    grid::Partitioner partitioner( p );
    if ( not partitioner ) {
        if ( config.has( "partitioner" ) ) {
            partitioner = grid::Partitioner( config.getSubConfiguration( "partitioner" ) );
        }
        else {
            if ( grid_->domain().global() ) {
                partitioner = grid::Partitioner( "equal_regions" );
            }
            else {
                partitioner = grid::Partitioner( "checkerboard" );
            }
        }
    }

    grid::Distribution distribution;
    ATLAS_TRACE_SCOPE( "Partitioning grid" ) { distribution = grid::Distribution( grid, partitioner ); }

    setup( distribution, config );
}

// ----------------------------------------------------------------------------

void StructuredColumns::compute_xy( idx_t i, idx_t j, PointXY& xy ) const {
    idx_t jj;
    if ( j < 0 ) {
        jj     = -j - 1 + north_pole_included_;
        xy[YY] = 180. - grid_->y( jj );
    }
    else if ( j >= ny_ ) {
        jj     = 2 * ny_ - j - 1 - south_pole_included_;
        xy[YY] = -180. - grid_->y( jj );
    }
    else {
        jj     = j;
        xy[YY] = grid_->y( jj );
    }
    xy[XX] = grid_->x( i, jj );
}


void StructuredColumns::Map2to1::resize( std::array<idx_t, 2> i_range, std::array<idx_t, 2> j_range ) {
    i_min_    = i_range[0];
    i_max_    = i_range[1];
    j_min_    = j_range[0];
    j_max_    = j_range[1];
    j_stride_ = ( i_max_ - i_min_ + 1 );
    data_.resize( ( i_max_ - i_min_ + 1 ) * ( j_max_ - j_min_ + 1 ) );
    atlas::omp::fill( data_.begin(), data_.end(), missing() + 1 );
}

const util::PartitionPolygon& StructuredColumns::polygon( idx_t halo ) const {
    if ( not polygon_ ) {
        polygon_.reset( new grid::StructuredPartitionPolygon( *this, halo ) );
    }
    return *polygon_;
}

const util::PartitionPolygons& StructuredColumns::polygons() const {
    if ( polygons_.size() ) {
        return polygons_;
    }
    polygon().allGather( polygons_ );
    return polygons_;
}

// ----------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------
StructuredColumns::~StructuredColumns() {
    delete grid_;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
Field StructuredColumns::createField( const eckit::Configuration& options ) const {
    Field field( config_name( options ), config_datatype( options ), config_shape( options ) );
    set_field_metadata( options, field );
    return field;
}

Field StructuredColumns::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | option::levels( other.levels() ) |
                        option::variables( other.variables() ) |
                        option::type( other.metadata().getString( "type", "scalar" ) ) | config );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const {
    ATLAS_ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& loc      = local_fieldset[f];
        Field& glb            = global_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> loc_field( make_leveled_view<const int>( loc ) );
            parallel::Field<int> glb_field( make_leveled_view<int>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> loc_field( make_leveled_view<const long>( loc ) );
            parallel::Field<long> glb_field( make_leveled_view<long>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> loc_field( make_leveled_view<const float>( loc ) );
            parallel::Field<float> glb_field( make_leveled_view<float>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> loc_field( make_leveled_view<const double>( loc ) );
            parallel::Field<double> glb_field( make_leveled_view<double>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else {
            throw_Exception( "datatype not supported", Here() );
        }
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void StructuredColumns::gather( const Field& local, Field& global ) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add( local );
    global_fields.add( global );
    gather( local_fields, global_fields );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const {
    ATLAS_ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& glb      = global_fieldset[f];
        Field& loc            = local_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> glb_field( make_leveled_view<const int>( glb ) );
            parallel::Field<int> loc_field( make_leveled_view<int>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> glb_field( make_leveled_view<const long>( glb ) );
            parallel::Field<long> loc_field( make_leveled_view<long>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> glb_field( make_leveled_view<const float>( glb ) );
            parallel::Field<float> loc_field( make_leveled_view<float>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> glb_field( make_leveled_view<const double>( glb ) );
            parallel::Field<double> loc_field( make_leveled_view<double>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else {
            throw_Exception( "datatype not supported", Here() );
        }

        glb.metadata().broadcast( loc.metadata(), root );
        loc.metadata().set( "global", false );
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void StructuredColumns::scatter( const Field& global, Field& local ) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add( global );
    local_fields.add( local );
    scatter( global_fields, local_fields );
}
// ----------------------------------------------------------------------------

std::string StructuredColumns::checksum( const FieldSet& fieldset ) const {
    eckit::MD5 md5;
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        const Field& field = fieldset[f];
        if ( field.datatype() == array::DataType::kind<int>() ) {
            md5 << checksum_3d_field<int>( checksum(), field );
        }
        else if ( field.datatype() == array::DataType::kind<long>() ) {
            md5 << checksum_3d_field<long>( checksum(), field );
        }
        else if ( field.datatype() == array::DataType::kind<float>() ) {
            md5 << checksum_3d_field<float>( checksum(), field );
        }
        else if ( field.datatype() == array::DataType::kind<double>() ) {
            md5 << checksum_3d_field<double>( checksum(), field );
        }
        else {
            throw_Exception( "datatype not supported", Here() );
        }
    }
    return md5;
}
std::string StructuredColumns::checksum( const Field& field ) const {
    FieldSet fieldset;
    fieldset.add( field );
    return checksum( fieldset );
}

const StructuredGrid& StructuredColumns::grid() const {
    return *grid_;
}

namespace {


template <int RANK>
struct FixupHaloForVectors {
    FixupHaloForVectors( const StructuredColumns& ) {}
    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector " ) {
            ATLAS_NOTIMPLEMENTED;
        }
    }
};

template <>
struct FixupHaloForVectors<2> {
    static constexpr int RANK = 2;
    const StructuredColumns& fs;
    FixupHaloForVectors( const StructuredColumns& _fs ) : fs( _fs ) {}

    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector" ) {
            auto array = array::make_view<DATATYPE, RANK>( field );
            for ( idx_t j = fs.j_begin_halo(); j < 0; ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n        = fs.index( i, j );
                    array( n, XX ) = -array( n, XX );
                    array( n, YY ) = -array( n, YY );
                }
            }
            for ( idx_t j = fs.grid().ny(); j < fs.j_end_halo(); ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n        = fs.index( i, j );
                    array( n, XX ) = -array( n, XX );
                    array( n, YY ) = -array( n, YY );
                }
            }
        }
    }
};

template <>
struct FixupHaloForVectors<3> {
    static constexpr int RANK = 3;
    const StructuredColumns& fs;
    FixupHaloForVectors( const StructuredColumns& _fs ) : fs( _fs ) {}

    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector" ) {
            auto array = array::make_view<DATATYPE, RANK>( field );
            for ( idx_t j = fs.j_begin_halo(); j < 0; ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n = fs.index( i, j );
                    for ( idx_t k = fs.k_begin(); k < fs.k_end(); ++k ) {
                        array( n, k, XX ) = -array( n, k, XX );
                        array( n, k, YY ) = -array( n, k, YY );
                    }
                }
            }
            for ( idx_t j = fs.grid().ny(); j < fs.j_end_halo(); ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n = fs.index( i, j );
                    for ( idx_t k = fs.k_begin(); k < fs.k_end(); ++k ) {
                        array( n, k, XX ) = -array( n, k, XX );
                        array( n, k, YY ) = -array( n, k, YY );
                    }
                }
            }
        }
    }
};


template <int RANK>
void dispatch_haloExchange( Field& field, const parallel::HaloExchange& halo_exchange, const StructuredColumns& fs ) {
    FixupHaloForVectors<RANK> fixup_halos( fs );
    if ( field.datatype() == array::DataType::kind<int>() ) {
        halo_exchange.template execute<int, RANK>( field.array(), false );
        fixup_halos.template apply<int>( field );
    }
    else if ( field.datatype() == array::DataType::kind<long>() ) {
        halo_exchange.template execute<long, RANK>( field.array(), false );
        fixup_halos.template apply<long>( field );
    }
    else if ( field.datatype() == array::DataType::kind<float>() ) {
        halo_exchange.template execute<float, RANK>( field.array(), false );
        fixup_halos.template apply<float>( field );
    }
    else if ( field.datatype() == array::DataType::kind<double>() ) {
        halo_exchange.template execute<double, RANK>( field.array(), false );
        fixup_halos.template apply<double>( field );
    }
    else {
        throw_Exception( "datatype not supported", Here() );
    }
    field.set_dirty( false );
}


template <int RANK>
void dispatch_adjointHaloExchange( Field& field, const parallel::HaloExchange& halo_exchange,
                                   const StructuredColumns& fs ) {
    FixupHaloForVectors<RANK> fixup_halos( fs );
    if ( field.datatype() == array::DataType::kind<int>() ) {
        halo_exchange.template execute_adjoint<int, RANK>( field.array(), false );
        fixup_halos.template apply<int>( field );
    }
    else if ( field.datatype() == array::DataType::kind<long>() ) {
        halo_exchange.template execute_adjoint<long, RANK>( field.array(), false );
        fixup_halos.template apply<long>( field );
    }
    else if ( field.datatype() == array::DataType::kind<float>() ) {
        halo_exchange.template execute_adjoint<float, RANK>( field.array(), false );
        fixup_halos.template apply<float>( field );
    }
    else if ( field.datatype() == array::DataType::kind<double>() ) {
        halo_exchange.template execute_adjoint<double, RANK>( field.array(), false );
        fixup_halos.template apply<double>( field );
    }
    else {
        throw_Exception( "datatype not supported", Here() );
    }
    field.set_dirty( false );
}
}  // namespace

void StructuredColumns::haloExchange( const FieldSet& fieldset, bool ) const {
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        Field& field = const_cast<FieldSet&>( fieldset )[f];
        switch ( field.rank() ) {
            case 1:
                dispatch_haloExchange<1>( field, halo_exchange(), *this );
                break;
            case 2:
                dispatch_haloExchange<2>( field, halo_exchange(), *this );
                break;
            case 3:
                dispatch_haloExchange<3>( field, halo_exchange(), *this );
                break;
            case 4:
                dispatch_haloExchange<4>( field, halo_exchange(), *this );
                break;
            default:
                throw_Exception( "Rank not supported", Here() );
        }
    }
}

void StructuredColumns::adjointHaloExchange( const FieldSet& fieldset, bool ) const {
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        Field& field = const_cast<FieldSet&>( fieldset )[f];
        switch ( field.rank() ) {
            case 1:
                dispatch_adjointHaloExchange<1>( field, halo_exchange(), *this );
                break;
            case 2:
                dispatch_adjointHaloExchange<2>( field, halo_exchange(), *this );
                break;
            case 3:
                dispatch_adjointHaloExchange<3>( field, halo_exchange(), *this );
                break;
            case 4:
                dispatch_adjointHaloExchange<4>( field, halo_exchange(), *this );
                break;
            default:
                throw_Exception( "Rank not supported", Here() );
        }
    }
}

void StructuredColumns::haloExchange( const Field& field, bool ) const {
    FieldSet fieldset;
    fieldset.add( field );
    haloExchange( fieldset );
}

void StructuredColumns::adjointHaloExchange( const Field& field, bool ) const {
    FieldSet fieldset;
    fieldset.add( field );
    adjointHaloExchange( fieldset );
}

size_t StructuredColumns::footprint() const {
    size_t size = sizeof( *this );
    size += ij2gp_.footprint();
    if ( field_xy_ ) {
        size += field_xy_.footprint();
    }
    if ( field_partition_ ) {
        size += field_partition_.footprint();
    }
    if ( field_global_index_ ) {
        size += field_global_index_.footprint();
    }
    if ( field_remote_index_ ) {
        size += field_remote_index_.footprint();
    }
    if ( field_index_i_ ) {
        size += field_index_i_.footprint();
    }
    if ( field_index_j_ ) {
        size += field_index_j_.footprint();
    }
    return size;
}

// ----------------------------------------------------------------------------


}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
