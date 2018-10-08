/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <functional>
#include <limits>

#include "eckit/utils/MD5.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/detail/Cache.h"

#undef atlas_omp_critical_ordered
#define atlas_omp_critical_ordered atlas_omp_critical

#if ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif

namespace atlas {
namespace functionspace {
namespace detail {

namespace {

template <typename T>
array::LocalView<T, 3> make_leveled_view( const Field& field ) {
    using namespace array;
    if ( field.levels() ) {
        if ( field.variables() ) { return make_view<T, 3>( field ).slice( Range::all(), Range::all(), Range::all() ); }
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
array::LocalView<T, 2> make_leveled_scalar_view( const Field& field ) {
    using namespace array;
    if ( field.levels() )
        return make_view<T, 2>( field ).slice( Range::all(), Range::all() );
    else
        return make_view<T, 1>( field ).slice( Range::all(), Range::dummy() );
}

template <typename T>
array::LocalView<T, 2> make_surface_view( const Field& field ) {
    using namespace array;
    if ( field.variables() )
        return make_view<T, 2>( field ).slice( Range::all(), Range::all() );
    else
        return make_view<T, 1>( field ).slice( Range::all(), Range::dummy() );
}

template <typename T>
array::LocalView<T, 2> make_per_level_view( const Field& field ) {
    using namespace array;
    if ( field.rank() == 2 )
        return make_view<T, 2>( field ).slice( Range::all(), Range::all() );
    else
        return make_view<T, 1>( field ).slice( Range::all(), Range::dummy() );
}

}  // namespace

class NodeColumnsHaloExchangeCache : public util::Cache<std::string, parallel::HaloExchange>,
                                     public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::HaloExchange>;
    NodeColumnsHaloExchangeCache() : Base( "NodeColumnsHaloExchangeCache" ) {}

public:
    static NodeColumnsHaloExchangeCache& instance() {
        static NodeColumnsHaloExchangeCache inst;
        return inst;
    }
    eckit::SharedPtr<value_type> get_or_create( const Mesh& mesh, long halo ) {
        creator_type creator = std::bind( &NodeColumnsHaloExchangeCache::create, mesh, halo );
        return Base::get_or_create( key( *mesh.get(), halo ), creator );
    }
    virtual void onMeshDestruction( mesh::detail::MeshImpl& mesh ) {
        for ( long jhalo = 1; jhalo < mesh::Halo( mesh ).size(); ++jhalo ) {
            remove( key( mesh, jhalo ) );
        }
    }

private:
    static Base::key_type key( const mesh::detail::MeshImpl& mesh, long halo ) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "],halo[size=" << halo << "]";
        return key.str();
    }

    static value_type* create( const Mesh& mesh, long halo ) {
        mesh.get()->attachObserver( instance() );

        value_type* value = new value_type();

        std::ostringstream ss;
        ss << "nb_nodes_including_halo[" << halo << "]";
        idx_t nb_nodes( mesh.nodes().size() );
        mesh.metadata().get( ss.str(), nb_nodes );

        value->setup( array::make_view<int, 1>( mesh.nodes().partition() ).data(),
                      array::make_view<idx_t, 1>( mesh.nodes().remote_index() ).data(), REMOTE_IDX_BASE, nb_nodes );

        return value;
    }
};

class NodeColumnsGatherScatterCache : public util::Cache<std::string, parallel::GatherScatter>,
                                      public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::GatherScatter>;
    NodeColumnsGatherScatterCache() : Base( "NodeColumnsGatherScatterCache" ) {}

public:
    static NodeColumnsGatherScatterCache& instance() {
        static NodeColumnsGatherScatterCache inst;
        return inst;
    }
    eckit::SharedPtr<value_type> get_or_create( const Mesh& mesh ) {
        creator_type creator = std::bind( &NodeColumnsGatherScatterCache::create, mesh );
        return Base::get_or_create( key( *mesh.get() ), creator );
    }
    virtual void onMeshDestruction( mesh::detail::MeshImpl& mesh ) { remove( key( mesh ) ); }

private:
    static Base::key_type key( const mesh::detail::MeshImpl& mesh ) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "]";
        return key.str();
    }

    static value_type* create( const Mesh& mesh ) {
        mesh.get()->attachObserver( instance() );

        value_type* value = new value_type();

        mesh::IsGhostNode is_ghost( mesh.nodes() );
        std::vector<int> mask( mesh.nodes().size() );
        const idx_t npts = mask.size();
        atlas_omp_parallel_for( idx_t n = 0; n < npts; ++n ) {
            mask[n] = is_ghost( n ) ? 1 : 0;

            // --> This would add periodic west-bc to the gather, but means that
            // global-sums, means, etc are computed wrong
            // if( mask[j] == 1 &&
            // internals::Topology::check(flags(j),internals::Topology::BC) ) {
            //  mask[j] = 0;
            //}
        }

        value->setup( array::make_view<int, 1>( mesh.nodes().partition() ).data(),
                      array::make_view<idx_t, 1>( mesh.nodes().remote_index() ).data(), REMOTE_IDX_BASE,
                      array::make_view<gidx_t, 1>( mesh.nodes().global_index() ).data(), mask.data(),
                      mesh.nodes().size() );
        return value;
    }
};

class NodeColumnsChecksumCache : public util::Cache<std::string, parallel::Checksum>,
                                 public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::Checksum>;
    NodeColumnsChecksumCache() : Base( "NodeColumnsChecksumCache" ) {}

public:
    static NodeColumnsChecksumCache& instance() {
        static NodeColumnsChecksumCache inst;
        return inst;
    }
    eckit::SharedPtr<value_type> get_or_create( const Mesh& mesh ) {
        creator_type creator = std::bind( &NodeColumnsChecksumCache::create, mesh );
        return Base::get_or_create( key( *mesh.get() ), creator );
    }
    virtual void onMeshDestruction( mesh::detail::MeshImpl& mesh ) { remove( key( mesh ) ); }

private:
    static Base::key_type key( const mesh::detail::MeshImpl& mesh ) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "]";
        return key.str();
    }

    static value_type* create( const Mesh& mesh ) {
        mesh.get()->attachObserver( instance() );
        value_type* value = new value_type();
        eckit::SharedPtr<parallel::GatherScatter> gather(
            NodeColumnsGatherScatterCache::instance().get_or_create( mesh ) );
        value->setup( gather );
        return value;
    }
};

NodeColumns::NodeColumns( Mesh mesh ) : NodeColumns( mesh, util::NoConfig() ) {}

NodeColumns::NodeColumns( Mesh mesh, const eckit::Configuration& config ) :
    mesh_( mesh ),
    nodes_( mesh_.nodes() ),
    nb_levels_( config.getInt( "levels", 0 ) ),
    nb_nodes_( 0 ) {
    ATLAS_TRACE();
    if ( config.has( "halo" ) ) { halo_ = mesh::Halo( config.getInt( "halo" ) ); }
    else {
        halo_ = mesh::Halo( mesh );
    }
    ASSERT( mesh_ );
    mesh::actions::build_nodes_parallel_fields( mesh_.nodes() );
    mesh::actions::build_periodic_boundaries( mesh_ );

    if ( halo_.size() > 0 ) {
        mesh::actions::build_halo( mesh_, halo_.size() );
        std::stringstream ss;
        ss << "nb_nodes_including_halo[" << halo_.size() << "]";
        mesh_.metadata().get( ss.str(), nb_nodes_ );
    }
    if ( !nb_nodes_ ) {
        std::stringstream ss;
        ss << "nb_nodes_including_halo[" << halo_.size() << "]";
        if ( !mesh_.metadata().get( ss.str(), nb_nodes_ ) ) { nb_nodes_ = mesh_.nodes().size(); }
    }
}

NodeColumns::~NodeColumns() {}

std::string NodeColumns::distribution() const {
    return mesh().metadata().getString( "distribution" );
}

idx_t NodeColumns::nb_nodes() const {
    return nb_nodes_;
}

idx_t NodeColumns::nb_nodes_global() const {
    if ( nb_nodes_global_ >= 0 ) return nb_nodes_global_;
    if ( Grid grid = mesh().grid() ) { nb_nodes_global_ = grid.size(); }
    else {
        nb_nodes_global_ = gather().glb_dof();
    }
    return nb_nodes_global_;
}

idx_t NodeColumns::config_nb_nodes( const eckit::Configuration& config ) const {
    idx_t size = nb_nodes();
    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            idx_t _nb_nodes_global = nb_nodes_global();
            size                   = ( mpi::comm().rank() == owner ? _nb_nodes_global : 0 );
        }
    }
    return size;
}

void NodeColumns::set_field_metadata( const eckit::Configuration& config, Field& field ) const {
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
}

array::DataType NodeColumns::config_datatype( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) throw eckit::AssertionFailed( "datatype missing", Here() );
    return array::DataType( kind );
}

std::string NodeColumns::config_name( const eckit::Configuration& config ) const {
    std::string name;
    config.get( "name", name );
    return name;
}

idx_t NodeColumns::config_levels( const eckit::Configuration& config ) const {
    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    return levels;
}

array::ArrayShape NodeColumns::config_shape( const eckit::Configuration& config ) const {
    array::ArrayShape shape;

    shape.push_back( config_nb_nodes( config ) );

    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    if ( levels > 0 ) shape.push_back( levels );

    idx_t variables( 0 );
    config.get( "variables", variables );
    if ( variables > 0 ) shape.push_back( variables );

    return shape;
}

Field NodeColumns::createField( const eckit::Configuration& config ) const {
    Field field = Field( config_name( config ), config_datatype( config ), config_shape( config ) );

    set_field_metadata( config, field );

    return field;
}

Field NodeColumns::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | option::levels( other.levels() ) |
                        option::variables( other.variables() ) | config );
}

namespace {

template <typename DATA_TYPE, int RANK>
array::ArrayView<DATA_TYPE, RANK> get_field_view( const Field& field, bool on_device ) {
    return on_device ? array::make_device_view<DATA_TYPE, RANK>( field ) : array::make_view<DATA_TYPE, RANK>( field );
}

template <int RANK>
void dispatch_haloExchange( Field& field, const parallel::HaloExchange& halo_exchange, bool on_device ) {
    if ( field.datatype() == array::DataType::kind<int>() ) {
        halo_exchange.template execute<int, RANK>( field.array(), on_device );
    }
    else if ( field.datatype() == array::DataType::kind<long>() ) {
        halo_exchange.template execute<long, RANK>( field.array(), on_device );
    }
    else if ( field.datatype() == array::DataType::kind<float>() ) {
        halo_exchange.template execute<float, RANK>( field.array(), on_device );
    }
    else if ( field.datatype() == array::DataType::kind<double>() ) {
        halo_exchange.template execute<double, RANK>( field.array(), on_device );
    }
    else
        throw eckit::Exception( "datatype not supported", Here() );
}
}  // namespace

void NodeColumns::haloExchange( FieldSet& fieldset, bool on_device ) const {
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        Field& field = fieldset[f];
        switch ( field.rank() ) {
            case 1:
                dispatch_haloExchange<1>( field, halo_exchange(), on_device );
                break;
            case 2:
                dispatch_haloExchange<2>( field, halo_exchange(), on_device );
                break;
            case 3:
                dispatch_haloExchange<3>( field, halo_exchange(), on_device );
                break;
            case 4:
                dispatch_haloExchange<4>( field, halo_exchange(), on_device );
                break;
            default:
                throw eckit::Exception( "Rank not supported", Here() );
        }
    }
}

void NodeColumns::haloExchange( Field& field, bool on_device ) const {
    FieldSet fieldset;
    fieldset.add( field );
    haloExchange( fieldset, on_device );
}
const parallel::HaloExchange& NodeColumns::halo_exchange() const {
    if ( halo_exchange_ ) return *halo_exchange_;
    halo_exchange_ = NodeColumnsHaloExchangeCache::instance().get_or_create( mesh_, halo_.size() );
    return *halo_exchange_;
}

void NodeColumns::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const {
    ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& loc      = local_fieldset[f];
        Field& glb            = global_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> loc_field( make_leveled_view<int>( loc ) );
            parallel::Field<int> glb_field( make_leveled_view<int>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> loc_field( make_leveled_view<long>( loc ) );
            parallel::Field<long> glb_field( make_leveled_view<long>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> loc_field( make_leveled_view<float>( loc ) );
            parallel::Field<float> glb_field( make_leveled_view<float>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> loc_field( make_leveled_view<double>( loc ) );
            parallel::Field<double> glb_field( make_leveled_view<double>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

void NodeColumns::gather( const Field& local, Field& global ) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add( local );
    global_fields.add( global );
    gather( local_fields, global_fields );
}
const parallel::GatherScatter& NodeColumns::gather() const {
    if ( gather_scatter_ ) return *gather_scatter_;
    gather_scatter_ = NodeColumnsGatherScatterCache::instance().get_or_create( mesh_ );
    return *gather_scatter_;
}
const parallel::GatherScatter& NodeColumns::scatter() const {
    if ( gather_scatter_ ) return *gather_scatter_;
    gather_scatter_ = NodeColumnsGatherScatterCache::instance().get_or_create( mesh_ );
    return *gather_scatter_;
}

void NodeColumns::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const {
    ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& glb      = global_fieldset[f];
        Field& loc            = local_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> glb_field( make_leveled_view<int>( glb ) );
            parallel::Field<int> loc_field( make_leveled_view<int>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> glb_field( make_leveled_view<long>( glb ) );
            parallel::Field<long> loc_field( make_leveled_view<long>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> glb_field( make_leveled_view<float>( glb ) );
            parallel::Field<float> loc_field( make_leveled_view<float>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> glb_field( make_leveled_view<double>( glb ) );
            parallel::Field<double> loc_field( make_leveled_view<double>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else
            throw eckit::Exception( "datatype not supported", Here() );

        glb.metadata().broadcast( loc.metadata(), root );
        loc.metadata().set( "global", false );
    }
}

void NodeColumns::scatter( const Field& global, Field& local ) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add( global );
    local_fields.add( local );
    scatter( global_fields, local_fields );
}

namespace {
template <typename T>
std::string checksum_3d_field( const parallel::Checksum& checksum, const Field& field ) {
    array::LocalView<T, 3> values = make_leveled_view<T>( field );
    array::ArrayT<T> surface_field( values.shape( 0 ), values.shape( 2 ) );
    array::ArrayView<T, 2> surface = array::make_view<T, 2>( surface_field );
    const idx_t npts               = values.shape( 0 );
    atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t j = 0; j < surface.shape( 1 ); ++j ) {
            surface( n, j ) = 0.;
            for ( idx_t l = 0; l < values.shape( 1 ); ++l )
                surface( n, j ) += values( n, l, j );
        }
    }
    return checksum.execute( surface.data(), surface_field.stride( 0 ) );
}
}  // namespace

std::string NodeColumns::checksum( const FieldSet& fieldset ) const {
    eckit::MD5 md5;
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        const Field& field = fieldset[f];
        if ( field.datatype() == array::DataType::kind<int>() )
            md5 << checksum_3d_field<int>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<long>() )
            md5 << checksum_3d_field<long>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<float>() )
            md5 << checksum_3d_field<float>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<double>() )
            md5 << checksum_3d_field<double>( checksum(), field );
        else
            throw eckit::Exception( "datatype not supported", Here() );
    }
    return md5;
}
std::string NodeColumns::checksum( const Field& field ) const {
    FieldSet fieldset;
    fieldset.add( field );
    return checksum( fieldset );
}

const parallel::Checksum& NodeColumns::checksum() const {
    if ( checksum_ ) return *checksum_;
    checksum_ = NodeColumnsChecksumCache::instance().get_or_create( mesh_ );
    return *checksum_;
}

// std::string NodesFunctionSpace::checksum( const FieldSet& fieldset ) const {
//  const parallel::Checksum& checksum = mesh_.checksum().get(checksum_name());

//  eckit::MD5 md5;
//  for( idx_t f=0; f<fieldset.size(); ++f ) {
//    const Field& field=fieldset[f];
//    if     ( field.datatype() == array::DataType::kind<int>() )
//      md5 << checksum.execute( field.data<int>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<long>() )
//      md5 << checksum.execute( field.data<long>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<float>() )
//      md5 << checksum.execute( field.data<float>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<double>() )
//      md5 << checksum.execute( field.data<double>(), field.stride(0) );
//    else throw eckit::Exception("datatype not supported",Here());
//  }
//  return md5;
//}
// std::string NodesFunctionSpace::checksum( const Field& field ) const {
//  FieldSet fieldset;
//  fieldset.add(field);
//  return checksum(fieldset);
//}

namespace {
inline double sqr( const double& val ) {
    return val * val;
}
}  // namespace

namespace detail {  // Collectives implementation

template <typename T>
void dispatch_sum( const NodeColumns& fs, const Field& field, T& result, idx_t& N ) {
    const mesh::IsGhostNode is_ghost( fs.nodes() );
    const array::LocalView<T, 2> arr = make_leveled_scalar_view<T>( field );
    T local_sum                      = 0;
    const idx_t npts                 = std::min<idx_t>( arr.shape( 0 ), fs.nb_nodes() );
  atlas_omp_pragma( omp parallel for default(shared) reduction(+:local_sum) )
  for( idx_t n=0; n<npts; ++n ) {
      if ( !is_ghost( n ) ) {
          for ( idx_t l = 0; l < arr.shape( 1 ); ++l )
              local_sum += arr( n, l );
      }
  }
  ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( local_sum, result, eckit::mpi::sum() ); }

  N = fs.nb_nodes_global() * arr.shape( 1 );
}

template <typename T>
void sum( const NodeColumns& fs, const Field& field, T& result, idx_t& N ) {
    if ( field.datatype() == array::DataType::kind<T>() ) { return dispatch_sum( fs, field, result, N ); }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                int tmp;
                dispatch_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_INT64: {
                long tmp;
                dispatch_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_REAL32: {
                float tmp;
                dispatch_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_REAL64: {
                double tmp;
                dispatch_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void dispatch_sum( const NodeColumns& fs, const Field& field, std::vector<T>& result, idx_t& N ) {
    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    const mesh::IsGhostNode is_ghost( fs.nodes() );
    const idx_t nvar = arr.shape( 2 );
    std::vector<T> local_sum( nvar, 0 );
    result.resize( nvar );

    atlas_omp_parallel {
        std::vector<T> local_sum_private( nvar, 0 );
        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            if ( !is_ghost( n ) ) {
                for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                    for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                        local_sum_private[j] += arr( n, l, j );
                    }
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t j = 0; j < nvar; ++j ) {
                local_sum[j] += local_sum_private[j];
            }
        }
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( local_sum, result, eckit::mpi::sum() ); }

    N = fs.nb_nodes_global() * arr.shape( 1 );
}

template <typename T>
void sum( const NodeColumns& fs, const Field& field, std::vector<T>& result, idx_t& N ) {
    if ( field.datatype() == array::DataType::kind<T>() ) { return dispatch_sum( fs, field, result, N ); }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void dispatch_sum_per_level( const NodeColumns& fs, const Field& field, Field& sum, idx_t& N ) {
    mesh::IsGhostNode is_ghost( fs.nodes() );

    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    sum.resize( shape );

    auto arr = make_leveled_view<T>( field );

    auto sum_per_level = make_per_level_view<T>( sum );

    for ( idx_t l = 0; l < sum_per_level.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < sum_per_level.shape( 1 ); ++j ) {
            sum_per_level( l, j ) = 0;
        }
    }

    atlas_omp_parallel {
        array::ArrayT<T> sum_per_level_private( sum_per_level.shape( 0 ), sum_per_level.shape( 1 ) );
        array::ArrayView<T, 2> sum_per_level_private_view = array::make_view<T, 2>( sum_per_level_private );

        for ( idx_t l = 0; l < sum_per_level_private_view.shape( 0 ); ++l ) {
            for ( idx_t j = 0; j < sum_per_level_private_view.shape( 1 ); ++j ) {
                sum_per_level_private_view( l, j ) = 0;
            }
        }

        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            if ( !is_ghost( n ) ) {
                for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                    for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                        sum_per_level_private_view( l, j ) += arr( n, l, j );
                    }
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t l = 0; l < sum_per_level_private.shape( 0 ); ++l ) {
                for ( idx_t j = 0; j < sum_per_level_private.shape( 1 ); ++j ) {
                    sum_per_level( l, j ) += sum_per_level_private_view( l, j );
                }
            }
        }
    }
    ATLAS_TRACE_MPI( ALLREDUCE ) {
        mpi::comm().allReduceInPlace( sum_per_level.data(), sum.size(), eckit::mpi::sum() );
    }
    N = fs.nb_nodes_global();
}

void sum_per_level( const NodeColumns& fs, const Field& field, Field& sum, idx_t& N ) {
    if ( field.datatype() != sum.datatype() ) {
        throw eckit::Exception( "Field and sum are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_sum_per_level<int>( fs, field, sum, N );
        case array::DataType::KIND_INT64:
            return dispatch_sum_per_level<long>( fs, field, sum, N );
        case array::DataType::KIND_REAL32:
            return dispatch_sum_per_level<float>( fs, field, sum, N );
        case array::DataType::KIND_REAL64:
            return dispatch_sum_per_level<double>( fs, field, sum, N );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename DATATYPE>
void dispatch_order_independent_sum_2d( const NodeColumns& fs, const Field& field, DATATYPE& result, idx_t& N ) {
    idx_t root   = 0;
    Field global = fs.createField( field, option::global() );
    fs.gather( field, global );
    result   = 0;
    auto glb = array::make_view<DATATYPE, 1>( global );
    for ( idx_t jnode = 0; jnode < glb.size(); ++jnode ) {
        result += glb( jnode );
    }
    ATLAS_TRACE_MPI( BROADCAST ) { mpi::comm().broadcast( &result, 1, root ); }
    N = fs.nb_nodes_global();
}

template <typename T>
void dispatch_order_independent_sum( const NodeColumns& fs, const Field& field, T& result, idx_t& N ) {
    if ( field.levels() ) {
        const array::LocalView<T, 2> arr = make_leveled_scalar_view<T>( field );

        Field surface_field = fs.createField<T>( option::name( "surface" ) | option::levels( false ) );
        auto surface        = array::make_view<T, 1>( surface_field );

        const idx_t N0 = std::min( surface.shape( 0 ), arr.shape( 0 ) );
        const idx_t N1 = arr.shape( 1 );
        for ( idx_t n = 0; n < N0; ++n ) {
            surface( n ) = 0;
            for ( idx_t l = 0; l < N1; ++l ) {
                surface( n ) += arr( n, l );
            }
        }
        dispatch_order_independent_sum_2d( fs, surface_field, result, N );
        N *= arr.shape( 1 );
    }
    else {
        dispatch_order_independent_sum_2d( fs, field, result, N );
    }
}

template <typename T>
void order_independent_sum( const NodeColumns& fs, const Field& field, T& result, idx_t& N ) {
    if ( field.datatype() == array::DataType::kind<T>() ) {
        return dispatch_order_independent_sum( fs, field, result, N );
    }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                int tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_INT64: {
                long tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_REAL32: {
                float tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            case array::DataType::KIND_REAL64: {
                double tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result = tmp;
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename DATATYPE>
void dispatch_order_independent_sum_2d( const NodeColumns& fs, const Field& field, std::vector<DATATYPE>& result,
                                        idx_t& N ) {
    idx_t nvar = field.variables();
    result.resize( nvar );
    for ( idx_t j = 0; j < nvar; ++j )
        result[j] = 0.;
    Field global = fs.createField( field, option::name( "global" ) | option::global() );
    fs.gather( field, global );
    if ( mpi::comm().rank() == 0 ) {
        const auto glb = make_surface_view<DATATYPE>( global );
        for ( idx_t n = 0; n < fs.nb_nodes_global(); ++n ) {
            for ( idx_t j = 0; j < nvar; ++j ) {
                result[j] += glb( n, j );
            }
        }
    }
    idx_t root = global.metadata().get<idx_t>( "owner" );
    ATLAS_TRACE_MPI( BROADCAST ) { mpi::comm().broadcast( result, root ); }
    N = fs.nb_nodes_global();
}

template <typename T>
void dispatch_order_independent_sum( const NodeColumns& fs, const Field& field, std::vector<T>& result, idx_t& N ) {
    if ( field.levels() ) {
        const idx_t nvar = field.variables();
        const auto arr   = make_leveled_view<T>( field );

        Field surface_field =
            fs.createField<T>( option::name( "surface" ) | option::variables( nvar ) | option::levels( false ) );
        auto surface = make_surface_view<T>( surface_field );

        atlas_omp_for( idx_t n = 0; n < arr.shape( 0 ); ++n ) {
            for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                surface( n, j ) = 0;
            }
        }

        for ( idx_t n = 0; n < arr.shape( 0 ); ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    surface( n, j ) += arr( n, l, j );
                }
            }
        }

        dispatch_order_independent_sum_2d( fs, surface_field, result, N );
        N *= arr.shape( 1 );
    }
    else {
        dispatch_order_independent_sum_2d( fs, field, result, N );
    }
}

template <typename T>
void order_independent_sum( const NodeColumns& fs, const Field& field, std::vector<T>& result, idx_t& N ) {
    if ( field.datatype() == array::DataType::kind<T>() ) {
        return dispatch_order_independent_sum( fs, field, result, N );
    }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_order_independent_sum( fs, field, tmp, N );
                result.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void dispatch_order_independent_sum_per_level( const NodeColumns& fs, const Field& field, Field& sumfield, idx_t& N ) {
    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    sumfield.resize( shape );

    auto sum = make_per_level_view<T>( sumfield );
    for ( idx_t l = 0; l < sum.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < sum.shape( 1 ); ++j ) {
            sum( l, j ) = 0.;
        }
    }

    idx_t root   = 0;
    Field global = fs.createField( field, option::name( "global" ) | option::global() );

    fs.gather( field, global );
    if ( mpi::comm().rank() == 0 ) {
        const array::LocalView<T, 3> glb = make_leveled_view<T>( global );

        for ( idx_t n = 0; n < glb.shape( 0 ); ++n ) {
            for ( idx_t l = 0; l < glb.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < glb.shape( 2 ); ++j ) {
                    sum( l, j ) += glb( n, l, j );
                }
            }
        }
    }
    ATLAS_TRACE_MPI( BROADCAST ) {
        std::vector<T> sum_array( sumfield.size() );
        if ( mpi::comm().rank() == root ) {
            idx_t c( 0 );
            for ( idx_t l = 0; l < sum.shape( 0 ); ++l ) {
                for ( idx_t j = 0; j < sum.shape( 1 ); ++j ) {
                    sum_array[c++] = sum( l, j );
                }
            }
        }
        mpi::comm().broadcast( sum_array, root );
        if ( mpi::comm().rank() != root ) {
            idx_t c( 0 );
            for ( idx_t l = 0; l < sum.shape( 0 ); ++l ) {
                for ( idx_t j = 0; j < sum.shape( 1 ); ++j ) {
                    sum( l, j ) = sum_array[c++];
                }
            }
        }
    }
    N = fs.nb_nodes_global();
}

void order_independent_sum_per_level( const NodeColumns& fs, const Field& field, Field& sum, idx_t& N ) {
    if ( field.datatype() != sum.datatype() ) {
        throw eckit::Exception( "Field and sum are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_order_independent_sum_per_level<int>( fs, field, sum, N );
        case array::DataType::KIND_INT64:
            return dispatch_order_independent_sum_per_level<long>( fs, field, sum, N );
        case array::DataType::KIND_REAL32:
            return dispatch_order_independent_sum_per_level<float>( fs, field, sum, N );
        case array::DataType::KIND_REAL64:
            return dispatch_order_independent_sum_per_level<double>( fs, field, sum, N );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void dispatch_minimum( const NodeColumns& fs, const Field& field, std::vector<T>& min ) {
    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    const idx_t nvar                 = arr.shape( 2 );
    min.resize( nvar );
    std::vector<T> local_minimum( nvar, std::numeric_limits<T>::max() );
    atlas_omp_parallel {
        std::vector<T> local_minimum_private( nvar, std::numeric_limits<T>::max() );
        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    local_minimum_private[j] = std::min( arr( n, l, j ), local_minimum_private[j] );
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                local_minimum[j] = std::min( local_minimum_private[j], local_minimum[j] );
            }
        }
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( local_minimum, min, eckit::mpi::min() ); }
}

template <typename T>
void minimum( const NodeColumns& fs, const Field& field, std::vector<T>& min ) {
    if ( field.datatype() == array::DataType::kind<T>() ) { return dispatch_minimum( fs, field, min ); }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_minimum( fs, field, tmp );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_minimum( fs, field, tmp );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_minimum( fs, field, tmp );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_minimum( fs, field, tmp );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void dispatch_maximum( const NodeColumns& fs, const Field& field, std::vector<T>& max ) {
    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    const idx_t nvar                 = arr.shape( 2 );
    max.resize( nvar );
    std::vector<T> local_maximum( nvar, -std::numeric_limits<T>::max() );
    atlas_omp_parallel {
        std::vector<T> local_maximum_private( nvar, -std::numeric_limits<T>::max() );
        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    local_maximum_private[j] = std::max( arr( n, l, j ), local_maximum_private[j] );
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t j = 0; j < nvar; ++j ) {
                local_maximum[j] = std::max( local_maximum_private[j], local_maximum[j] );
            }
        }
    }
    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( local_maximum, max, eckit::mpi::max() ); }
}

template <typename T>
void maximum( const NodeColumns& fs, const Field& field, std::vector<T>& max ) {
    if ( field.datatype() == array::DataType::kind<T>() ) { return dispatch_maximum( fs, field, max ); }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_maximum( fs, field, tmp );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_maximum( fs, field, tmp );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_maximum( fs, field, tmp );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_maximum( fs, field, tmp );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void minimum( const NodeColumns& fs, const Field& field, T& min ) {
    std::vector<T> v;
    minimum( fs, field, v );
    min = v[0];
}

template <typename T>
void maximum( const NodeColumns& fs, const Field& field, T& max ) {
    std::vector<T> v;
    maximum( fs, field, v );
    max = v[0];
}

template <typename T>
void dispatch_minimum_per_level( const NodeColumns& fs, const Field& field, Field& min_field ) {
    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    min_field.resize( shape );
    auto min = make_per_level_view<T>( min_field );

    for ( idx_t l = 0; l < min.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < min.shape( 1 ); ++j ) {
            min( l, j ) = std::numeric_limits<T>::max();
        }
    }

    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    atlas_omp_parallel {
        array::ArrayT<T> min_private( min.shape( 0 ), min.shape( 1 ) );
        array::ArrayView<T, 2> min_private_view = array::make_view<T, 2>( min_private );
        for ( idx_t l = 0; l < min.shape( 0 ); ++l ) {
            for ( idx_t j = 0; j < min.shape( 1 ); ++j ) {
                min_private_view( l, j ) = std::numeric_limits<T>::max();
            }
        }

        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    min_private_view( l, j ) = std::min( arr( n, l, j ), min_private_view( l, j ) );
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    min( l, j ) = std::min( min_private_view( l, j ), min( l, j ) );
                }
            }
        }
    }
    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduceInPlace( min.data(), min_field.size(), eckit::mpi::min() ); }
}

void minimum_per_level( const NodeColumns& fs, const Field& field, Field& min ) {
    if ( field.datatype() != min.datatype() ) {
        throw eckit::Exception( "Field and min are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_minimum_per_level<int>( fs, field, min );
        case array::DataType::KIND_INT64:
            return dispatch_minimum_per_level<long>( fs, field, min );
        case array::DataType::KIND_REAL32:
            return dispatch_minimum_per_level<float>( fs, field, min );
        case array::DataType::KIND_REAL64:
            return dispatch_minimum_per_level<double>( fs, field, min );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void dispatch_maximum_per_level( const NodeColumns& fs, const Field& field, Field& max_field ) {
    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    max_field.resize( shape );
    auto max = make_per_level_view<T>( max_field );

    for ( idx_t l = 0; l < max.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < max.shape( 1 ); ++j ) {
            max( l, j ) = -std::numeric_limits<T>::max();
        }
    }

    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    atlas_omp_parallel {
        array::ArrayT<T> max_private( max.shape( 0 ), max.shape( 1 ) );
        array::ArrayView<T, 2> max_private_view = array::make_view<T, 2>( max_private );

        for ( idx_t l = 0; l < max_private_view.shape( 0 ); ++l ) {
            for ( idx_t j = 0; j < max_private_view.shape( 1 ); ++j ) {
                max_private_view( l, j ) = -std::numeric_limits<T>::max();
            }
        }

        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    max_private_view( l, j ) = std::max( arr( n, l, j ), max_private_view( l, j ) );
                }
            }
        }
        atlas_omp_critical {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < arr.shape( 2 ); ++j ) {
                    max( l, j ) = std::max( max_private_view( l, j ), max( l, j ) );
                }
            }
        }
    }
    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduceInPlace( max.data(), max_field.size(), eckit::mpi::max() ); }
}

void maximum_per_level( const NodeColumns& fs, const Field& field, Field& max ) {
    if ( field.datatype() != max.datatype() ) {
        throw eckit::Exception( "Field and max are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_maximum_per_level<int>( fs, field, max );
        case array::DataType::KIND_INT64:
            return dispatch_maximum_per_level<long>( fs, field, max );
        case array::DataType::KIND_REAL32:
            return dispatch_maximum_per_level<float>( fs, field, max );
        case array::DataType::KIND_REAL64:
            return dispatch_maximum_per_level<double>( fs, field, max );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void dispatch_minimum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& min,
                                    std::vector<gidx_t>& glb_idx, std::vector<idx_t>& level ) {
    array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    idx_t nvar                 = arr.shape( 2 );
    min.resize( nvar );
    glb_idx.resize( nvar );
    level.resize( nvar );
    std::vector<T> local_minimum( nvar, std::numeric_limits<T>::max() );
    std::vector<idx_t> loc_node( nvar );
    std::vector<idx_t> loc_level( nvar );
    atlas_omp_parallel {
        std::vector<T> local_minimum_private( nvar, std::numeric_limits<T>::max() );
        std::vector<idx_t> loc_node_private( nvar );
        std::vector<idx_t> loc_level_private( nvar );
        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( arr( n, l, j ) < local_minimum_private[j] ) {
                        local_minimum_private[j] = arr( n, l, j );
                        loc_node_private[j]      = n;
                        loc_level_private[j]     = l;
                    }
                }
            }
        }
        atlas_omp_critical_ordered {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( local_minimum_private[j] < local_minimum[j] ) {
                        local_minimum[j] = local_minimum_private[j];
                        loc_node[j]      = loc_node_private[j];
                        loc_level[j]     = loc_level_private[j];
                    }
                }
            }
        }
    }
    std::vector<std::pair<T, int>> min_and_gidx_loc( nvar );
    std::vector<std::pair<T, int>> min_and_level_loc( nvar );
    std::vector<std::pair<T, int>> min_and_gidx_glb( nvar );
    std::vector<std::pair<T, int>> min_and_level_glb( nvar );
    const array::ArrayView<gidx_t, 1> global_index = array::make_view<gidx_t, 1>( fs.nodes().global_index() );
    for ( idx_t j = 0; j < nvar; ++j ) {
        gidx_t glb_idx = global_index( loc_node[j] );
        ASSERT( glb_idx < std::numeric_limits<int>::max() );  // pairs with 64bit
                                                              // integer for second not
                                                              // implemented
        min_and_gidx_loc[j]  = std::make_pair( local_minimum[j], glb_idx );
        min_and_level_loc[j] = std::make_pair( local_minimum[j], loc_level[j] );
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) {
        mpi::comm().allReduce( min_and_gidx_loc, min_and_gidx_glb, eckit::mpi::minloc() );
        mpi::comm().allReduce( min_and_level_loc, min_and_level_glb, eckit::mpi::minloc() );
    }

    for ( idx_t j = 0; j < nvar; ++j ) {
        min[j]     = min_and_gidx_glb[j].first;
        glb_idx[j] = min_and_gidx_glb[j].second;
        level[j]   = min_and_level_glb[j].second;
    }
}

template <typename T>
void minimum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx,
                           std::vector<idx_t>& level ) {
    if ( field.datatype() == array::DataType::kind<T>() ) {
        return dispatch_minimum_and_location( fs, field, min, glb_idx, level );
    }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_minimum_and_location( fs, field, tmp, glb_idx, level );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_minimum_and_location( fs, field, tmp, glb_idx, level );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_minimum_and_location( fs, field, tmp, glb_idx, level );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_minimum_and_location( fs, field, tmp, glb_idx, level );
                min.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void dispatch_maximum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& max,
                                    std::vector<gidx_t>& glb_idx, std::vector<idx_t>& level ) {
    array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    idx_t nvar                 = arr.shape( 2 );
    max.resize( nvar );
    glb_idx.resize( nvar );
    level.resize( nvar );
    std::vector<T> local_maximum( nvar, -std::numeric_limits<T>::max() );
    std::vector<idx_t> loc_node( nvar );
    std::vector<idx_t> loc_level( nvar );
    atlas_omp_parallel {
        std::vector<T> local_maximum_private( nvar, -std::numeric_limits<T>::max() );
        std::vector<idx_t> loc_node_private( nvar );
        std::vector<idx_t> loc_level_private( nvar );
        const idx_t npts = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( arr( n, l, j ) > local_maximum_private[j] ) {
                        local_maximum_private[j] = arr( n, l, j );
                        loc_node_private[j]      = n;
                        loc_level_private[j]     = l;
                    }
                }
            }
        }
        atlas_omp_critical_ordered {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( local_maximum_private[j] > local_maximum[j] ) {
                        local_maximum[j] = local_maximum_private[j];
                        loc_node[j]      = loc_node_private[j];
                        loc_level[j]     = loc_level_private[j];
                    }
                }
            }
        }
    }
    std::vector<std::pair<T, int>> max_and_gidx_loc( nvar );
    std::vector<std::pair<T, int>> max_and_level_loc( nvar );
    std::vector<std::pair<T, int>> max_and_gidx_glb( nvar );
    std::vector<std::pair<T, int>> max_and_level_glb( nvar );
    const array::ArrayView<gidx_t, 1> global_index = array::make_view<gidx_t, 1>( fs.nodes().global_index() );
    for ( idx_t j = 0; j < nvar; ++j ) {
        gidx_t glb_idx = global_index( loc_node[j] );
        ASSERT( glb_idx < std::numeric_limits<int>::max() );  // pairs with 64bit
                                                              // integer for second not
                                                              // implemented
        max_and_gidx_loc[j]  = std::make_pair( local_maximum[j], glb_idx );
        max_and_level_loc[j] = std::make_pair( local_maximum[j], loc_level[j] );
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) {
        mpi::comm().allReduce( max_and_gidx_loc, max_and_gidx_glb, eckit::mpi::maxloc() );
        mpi::comm().allReduce( max_and_level_loc, max_and_level_glb, eckit::mpi::maxloc() );
    }

    for ( idx_t j = 0; j < nvar; ++j ) {
        max[j]     = max_and_gidx_glb[j].first;
        glb_idx[j] = max_and_gidx_glb[j].second;
        level[j]   = max_and_level_glb[j].second;
    }
}

template <typename T>
void maximum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx,
                           std::vector<idx_t>& level ) {
    if ( field.datatype() == array::DataType::kind<T>() ) {
        return dispatch_maximum_and_location( fs, field, max, glb_idx, level );
    }
    else {
        switch ( field.datatype().kind() ) {
            case array::DataType::KIND_INT32: {
                std::vector<int> tmp;
                dispatch_maximum_and_location( fs, field, tmp, glb_idx, level );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_INT64: {
                std::vector<long> tmp;
                dispatch_maximum_and_location( fs, field, tmp, glb_idx, level );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL32: {
                std::vector<float> tmp;
                dispatch_maximum_and_location( fs, field, tmp, glb_idx, level );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            case array::DataType::KIND_REAL64: {
                std::vector<double> tmp;
                dispatch_maximum_and_location( fs, field, tmp, glb_idx, level );
                max.assign( tmp.begin(), tmp.end() );
                return;
            }
            default:
                throw eckit::Exception( "datatype not supported", Here() );
        }
    }
}

template <typename T>
void minimum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& min,
                           std::vector<gidx_t>& glb_idx ) {
    std::vector<idx_t> level;
    minimum_and_location( fs, field, min, glb_idx, level );
}

template <typename T>
void maximum_and_location( const NodeColumns& fs, const Field& field, std::vector<T>& max,
                           std::vector<gidx_t>& glb_idx ) {
    std::vector<idx_t> level;
    maximum_and_location( fs, field, max, glb_idx, level );
}

template <typename T>
void minimum_and_location( const NodeColumns& fs, const Field& field, T& min, gidx_t& glb_idx, idx_t& level ) {
    std::vector<T> minv;
    std::vector<gidx_t> gidxv;
    std::vector<idx_t> levelv;
    minimum_and_location( fs, field, minv, gidxv, levelv );
    min     = minv[0];
    glb_idx = gidxv[0];
    level   = levelv[0];
}

template <typename T>
void maximum_and_location( const NodeColumns& fs, const Field& field, T& max, gidx_t& glb_idx, idx_t& level ) {
    std::vector<T> maxv;
    std::vector<gidx_t> gidxv;
    std::vector<idx_t> levelv;
    maximum_and_location( fs, field, maxv, gidxv, levelv );
    max     = maxv[0];
    glb_idx = gidxv[0];
    level   = levelv[0];
}

template <typename T>
void minimum_and_location( const NodeColumns& fs, const Field& field, T& min, gidx_t& glb_idx ) {
    idx_t level;
    minimum_and_location( fs, field, min, glb_idx, level );
}

template <typename T>
void maximum_and_location( const NodeColumns& fs, const Field& field, T& max, gidx_t& glb_idx ) {
    idx_t level;
    maximum_and_location( fs, field, max, glb_idx, level );
}

template <typename T>
void dispatch_minimum_and_location_per_level( const NodeColumns& fs, const Field& field, Field& min_field,
                                              Field& glb_idx_field ) {
    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    min_field.resize( shape );
    glb_idx_field.resize( shape );
    const idx_t nvar = arr.shape( 2 );
    auto min         = make_per_level_view<T>( min_field );
    auto glb_idx     = make_per_level_view<gidx_t>( glb_idx_field );

    for ( idx_t l = 0; l < min.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < min.shape( 1 ); ++j ) {
            min( l, j ) = std::numeric_limits<T>::max();
        }
    }

    atlas_omp_parallel {
        array::ArrayT<T> min_private( min.shape( 0 ), min.shape( 1 ) );
        array::ArrayView<T, 2> min_private_view = array::make_view<T, 2>( min_private );

        for ( idx_t l = 0; l < min_private_view.shape( 0 ); ++l ) {
            for ( idx_t j = 0; j < min_private_view.shape( 1 ); ++j ) {
                min_private_view( l, j ) = std::numeric_limits<T>::max();
            }
        }

        array::ArrayT<gidx_t> glb_idx_private( glb_idx.shape( 0 ), glb_idx.shape( 1 ) );
        array::ArrayView<gidx_t, 2> glb_idx_private_view = array::make_view<gidx_t, 2>( glb_idx_private );
        const idx_t npts                                 = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( arr( n, l, j ) < min( l, j ) ) {
                        min_private_view( l, j )     = arr( n, l, j );
                        glb_idx_private_view( l, j ) = n;
                    }
                }
            }
        }
        atlas_omp_critical_ordered {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( min_private_view( l, j ) < min( l, j ) ) {
                        min( l, j )     = min_private_view( l, j );
                        glb_idx( l, j ) = glb_idx_private_view( l, j );
                    }
                }
            }
        }
    }
    const idx_t nlev = arr.shape( 1 );
    std::vector<std::pair<T, int>> min_and_gidx_loc( nlev * nvar );
    std::vector<std::pair<T, int>> min_and_gidx_glb( nlev * nvar );
    const array::ArrayView<gidx_t, 1> global_index = array::make_view<gidx_t, 1>( fs.nodes().global_index() );
    atlas_omp_parallel_for( idx_t l = 0; l < nlev; ++l ) {
        for ( idx_t j = 0; j < nvar; ++j ) {
            gidx_t gidx = global_index( glb_idx( l, j ) );
            ASSERT( gidx < std::numeric_limits<int>::max() );  // pairs with 64bit
                                                               // integer for second not
                                                               // implemented
            min_and_gidx_loc[j + nvar * l] = std::make_pair( min( l, j ), gidx );
        }
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( min_and_gidx_loc, min_and_gidx_glb, eckit::mpi::minloc() ); }

    atlas_omp_parallel_for( idx_t l = 0; l < nlev; ++l ) {
        for ( idx_t j = 0; j < nvar; ++j ) {
            min( l, j )     = min_and_gidx_glb[j + l * nvar].first;
            glb_idx( l, j ) = min_and_gidx_glb[j + l * nvar].second;
        }
    }
}

void minimum_and_location_per_level( const NodeColumns& fs, const Field& field, Field& min, Field& glb_idx ) {
    if ( field.datatype() != min.datatype() ) {
        throw eckit::Exception( "Field and min are not of same datatype.", Here() );
    }
    if ( glb_idx.datatype() != array::DataType::kind<gidx_t>() ) {
        throw eckit::Exception( "glb_idx Field is not of correct datatype", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_minimum_and_location_per_level<int>( fs, field, min, glb_idx );
        case array::DataType::KIND_INT64:
            return dispatch_minimum_and_location_per_level<long>( fs, field, min, glb_idx );
        case array::DataType::KIND_REAL32:
            return dispatch_minimum_and_location_per_level<float>( fs, field, min, glb_idx );
        case array::DataType::KIND_REAL64:
            return dispatch_minimum_and_location_per_level<double>( fs, field, min, glb_idx );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void dispatch_maximum_and_location_per_level( const NodeColumns& fs, const Field& field, Field& max_field,
                                              Field& glb_idx_field ) {
    const array::LocalView<T, 3> arr = make_leveled_view<T>( field );
    array::ArrayShape shape;
    shape.reserve( field.rank() - 1 );
    for ( idx_t j = 1; j < field.rank(); ++j )
        shape.push_back( field.shape( j ) );
    max_field.resize( shape );
    glb_idx_field.resize( shape );
    const idx_t nvar = arr.shape( 2 );
    auto max         = make_per_level_view<T>( max_field );
    auto glb_idx     = make_per_level_view<gidx_t>( glb_idx_field );

    for ( idx_t l = 0; l < max.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < max.shape( 1 ); ++j ) {
            max( l, j ) = -std::numeric_limits<T>::max();
        }
    }

    atlas_omp_parallel {
        array::ArrayT<T> max_private( max.shape( 0 ), max.shape( 1 ) );
        array::ArrayView<T, 2> max_private_view = array::make_view<T, 2>( max_private );

        for ( idx_t l = 0; l < max_private_view.shape( 0 ); ++l ) {
            for ( idx_t j = 0; j < max_private_view.shape( 1 ); ++j ) {
                max_private_view( l, j ) = -std::numeric_limits<T>::max();
            }
        }

        array::ArrayT<gidx_t> glb_idx_private( glb_idx.shape( 0 ), glb_idx.shape( 1 ) );
        array::ArrayView<gidx_t, 2> glb_idx_private_view = array::make_view<gidx_t, 2>( glb_idx_private );
        const idx_t npts                                 = arr.shape( 0 );
        atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( arr( n, l, j ) > max( l, j ) ) {
                        max_private_view( l, j )     = arr( n, l, j );
                        glb_idx_private_view( l, j ) = n;
                    }
                }
            }
        }
        atlas_omp_critical_ordered {
            for ( idx_t l = 0; l < arr.shape( 1 ); ++l ) {
                for ( idx_t j = 0; j < nvar; ++j ) {
                    if ( max_private_view( l, j ) > max( l, j ) ) {
                        max( l, j )     = max_private_view( l, j );
                        glb_idx( l, j ) = glb_idx_private_view( l, j );
                    }
                }
            }
        }
    }

    const idx_t nlev = arr.shape( 1 );
    std::vector<std::pair<T, int>> max_and_gidx_loc( nlev * nvar );
    std::vector<std::pair<T, int>> max_and_gidx_glb( nlev * nvar );
    const array::ArrayView<gidx_t, 1> global_index = array::make_view<gidx_t, 1>( fs.nodes().global_index() );
    atlas_omp_parallel_for( idx_t l = 0; l < nlev; ++l ) {
        for ( idx_t j = 0; j < nvar; ++j ) {
            gidx_t gidx = global_index( glb_idx( l, j ) );
            ASSERT( gidx < std::numeric_limits<int>::max() );  // pairs with 64bit
                                                               // integer for second not
                                                               // implemented
            max_and_gidx_loc[j + nvar * l] = std::make_pair( max( l, j ), gidx );
        }
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduce( max_and_gidx_loc, max_and_gidx_glb, eckit::mpi::maxloc() ); }

    atlas_omp_parallel_for( idx_t l = 0; l < nlev; ++l ) {
        for ( idx_t j = 0; j < nvar; ++j ) {
            max( l, j )     = max_and_gidx_glb[j + l * nvar].first;
            glb_idx( l, j ) = max_and_gidx_glb[j + l * nvar].second;
        }
    }
}

void maximum_and_location_per_level( const NodeColumns& fs, const Field& field, Field& max, Field& glb_idx ) {
    if ( field.datatype() != max.datatype() ) {
        throw eckit::Exception( "Field and max are not of same datatype.", Here() );
    }
    if ( glb_idx.datatype() != array::DataType::kind<gidx_t>() ) {
        throw eckit::Exception( "glb_idx Field is not of correct datatype", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_maximum_and_location_per_level<int>( fs, field, max, glb_idx );
        case array::DataType::KIND_INT64:
            return dispatch_maximum_and_location_per_level<long>( fs, field, max, glb_idx );
        case array::DataType::KIND_REAL32:
            return dispatch_maximum_and_location_per_level<float>( fs, field, max, glb_idx );
        case array::DataType::KIND_REAL64:
            return dispatch_maximum_and_location_per_level<double>( fs, field, max, glb_idx );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void mean( const NodeColumns& fs, const Field& field, T& result, idx_t& N ) {
    sum( fs, field, result, N );
    result /= static_cast<double>( N );
}

template <typename T>
void mean( const NodeColumns& fs, const Field& field, std::vector<T>& result, idx_t& N ) {
    sum( fs, field, result, N );
    for ( idx_t j = 0; j < result.size(); ++j ) {
        result[j] /= static_cast<double>( N );
    }
}

template <typename T>
void dispatch_mean_per_level( const NodeColumns& fs, const Field& field, Field& mean, idx_t& N ) {
    dispatch_sum_per_level<T>( fs, field, mean, N );
    auto view = make_per_level_view<T>( mean );
    for ( idx_t l = 0; l < view.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < view.shape( 1 ); ++j ) {
            view( l, j ) /= static_cast<double>( N );
        }
    }
}

void mean_per_level( const NodeColumns& fs, const Field& field, Field& mean, idx_t& N ) {
    if ( field.datatype() != mean.datatype() ) {
        throw eckit::Exception( "Field and sum are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_mean_per_level<int>( fs, field, mean, N );
        case array::DataType::KIND_INT64:
            return dispatch_mean_per_level<long>( fs, field, mean, N );
        case array::DataType::KIND_REAL32:
            return dispatch_mean_per_level<float>( fs, field, mean, N );
        case array::DataType::KIND_REAL64:
            return dispatch_mean_per_level<double>( fs, field, mean, N );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

template <typename T>
void mean_and_standard_deviation( const NodeColumns& fs, const Field& field, T& mu, T& sigma, idx_t& N ) {
    mean( fs, field, mu, N );
    Field squared_diff_field = fs.createField( option::name( "sqr_diff" ) | option::datatype( field.datatype() ) |
                                               option::levels( field.levels() ) );

    array::LocalView<T, 2> squared_diff = make_leveled_scalar_view<T>( squared_diff_field );
    array::LocalView<T, 2> values       = make_leveled_scalar_view<T>( field );

    const idx_t npts = std::min<idx_t>( values.shape( 0 ), fs.nb_nodes() );
    atlas_omp_parallel_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t l = 0; l < values.shape( 1 ); ++l ) {
            squared_diff( n, l ) = sqr( values( n, l ) - mu );
        }
    }
    mean( fs, squared_diff_field, sigma, N );
    sigma = std::sqrt( sigma );
}

template <typename T>
void mean_and_standard_deviation( const NodeColumns& fs, const Field& field, std::vector<T>& mu, std::vector<T>& sigma,
                                  idx_t& N ) {
    mean( fs, field, mu, N );
    Field squared_diff_field = fs.createField<T>( option::name( "sqr_diff" ) | option::levels( field.levels() ) |
                                                  option::variables( field.variables() ) );
    array::LocalView<T, 3> squared_diff = make_leveled_view<T>( squared_diff_field );
    array::LocalView<T, 3> values       = make_leveled_view<T>( field );

    const idx_t npts = std::min<idx_t>( values.shape( 0 ), fs.nb_nodes() );
    atlas_omp_parallel_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t l = 0; l < values.shape( 1 ); ++l ) {
            for ( idx_t j = 0; j < values.shape( 2 ); ++j ) {
                squared_diff( n, l, j ) = sqr( values( n, l, j ) - mu[j] );
            }
        }
    }
    mean( fs, squared_diff_field, sigma, N );
    for ( idx_t j = 0; j < sigma.size(); ++j ) {
        sigma[j] = std::sqrt( sigma[j] );
    }
}

template <typename T>
void dispatch_mean_and_standard_deviation_per_level( const NodeColumns& fs, const Field& field, Field& mean,
                                                     Field& stddev, idx_t& N ) {
    dispatch_mean_per_level<T>( fs, field, mean, N );
    Field squared_diff_field = fs.createField<T>( option::name( "sqr_diff" ) | option::levels( field.levels() ) |
                                                  option::variables( field.variables() ) );
    auto squared_diff        = make_leveled_view<T>( squared_diff_field );
    auto values              = make_leveled_view<T>( field );
    auto mu                  = make_per_level_view<T>( mean );

    const idx_t npts = std::min<idx_t>( values.shape( 0 ), fs.nb_nodes() );
    atlas_omp_parallel_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t l = 0; l < values.shape( 1 ); ++l ) {
            for ( idx_t j = 0; j < values.shape( 2 ); ++j ) {
                squared_diff( n, l, j ) = sqr( values( n, l, j ) - mu( l, j ) );
            }
        }
    }
    dispatch_mean_per_level<T>( fs, squared_diff_field, stddev, N );
    auto sigma = make_per_level_view<T>( stddev );
    atlas_omp_for( idx_t l = 0; l < sigma.shape( 0 ); ++l ) {
        for ( idx_t j = 0; j < sigma.shape( 1 ); ++j ) {
            sigma( l, j ) = std::sqrt( sigma( l, j ) );
        }
    }
}

void mean_and_standard_deviation_per_level( const NodeColumns& fs, const Field& field, Field& mean, Field& stddev,
                                            idx_t& N ) {
    if ( field.datatype() != mean.datatype() ) {
        throw eckit::Exception( "Field and mean are not of same datatype.", Here() );
    }
    if ( field.datatype() != stddev.datatype() ) {
        throw eckit::Exception( "Field and stddev are not of same datatype.", Here() );
    }
    switch ( field.datatype().kind() ) {
        case array::DataType::KIND_INT32:
            return dispatch_mean_and_standard_deviation_per_level<int>( fs, field, mean, stddev, N );
        case array::DataType::KIND_INT64:
            return dispatch_mean_and_standard_deviation_per_level<long>( fs, field, mean, stddev, N );
        case array::DataType::KIND_REAL32:
            return dispatch_mean_and_standard_deviation_per_level<float>( fs, field, mean, stddev, N );
        case array::DataType::KIND_REAL64:
            return dispatch_mean_and_standard_deviation_per_level<double>( fs, field, mean, stddev, N );
        default:
            throw eckit::Exception( "datatype not supported", Here() );
    }
}

}  // namespace detail

template <typename Value>
NodeColumns::FieldStatisticsT<Value>::FieldStatisticsT( const NodeColumns* f ) : functionspace( *f ) {}

template <typename Vector>
NodeColumns::FieldStatisticsVectorT<Vector>::FieldStatisticsVectorT( const NodeColumns* f ) : functionspace( *f ) {}

NodeColumns::FieldStatistics::FieldStatistics( const NodeColumns* f ) : functionspace( *f ) {}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::sum( const Field& field, Value& result, idx_t& N ) const {
    detail::sum( functionspace, field, result, N );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::sum( const Field& field, Vector& result, idx_t& N ) const {
    detail::sum( functionspace, field, result, N );
}

void NodeColumns::FieldStatistics::sumPerLevel( const Field& field, Field& result, idx_t& N ) const {
    detail::sum_per_level( functionspace, field, result, N );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::orderIndependentSum( const Field& field, Value& result, idx_t& N ) const {
    detail::order_independent_sum( functionspace, field, result, N );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::orderIndependentSum( const Field& field, Vector& result,
                                                                       idx_t& N ) const {
    detail::order_independent_sum( functionspace, field, result, N );
}

void NodeColumns::FieldStatistics::orderIndependentSumPerLevel( const Field& field, Field& result, idx_t& N ) const {
    detail::order_independent_sum_per_level( functionspace, field, result, N );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::minimum( const Field& field, Value& minimum ) const {
    detail::minimum( functionspace, field, minimum );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::maximum( const Field& field, Value& maximum ) const {
    detail::maximum( functionspace, field, maximum );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::minimum( const Field& field, Vector& minimum ) const {
    detail::minimum( functionspace, field, minimum );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::maximum( const Field& field, Vector& maximum ) const {
    detail::maximum( functionspace, field, maximum );
}

void NodeColumns::FieldStatistics::minimumPerLevel( const Field& field, Field& minimum ) const {
    detail::minimum_per_level( functionspace, field, minimum );
}

void NodeColumns::FieldStatistics::maximumPerLevel( const Field& field, Field& maximum ) const {
    detail::maximum_per_level( functionspace, field, maximum );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::minimumAndLocation( const Field& field, Value& minimum,
                                                               gidx_t& glb_idx ) const {
    detail::minimum_and_location( functionspace, field, minimum, glb_idx );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::maximumAndLocation( const Field& field, Value& maximum,
                                                               gidx_t& glb_idx ) const {
    detail::maximum_and_location( functionspace, field, maximum, glb_idx );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::minimumAndLocation( const Field& field, Value& minimum, gidx_t& glb_idx,
                                                               idx_t& level ) const {
    detail::minimum_and_location( functionspace, field, minimum, glb_idx, level );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::maximumAndLocation( const Field& field, Value& maximum, gidx_t& glb_idx,
                                                               idx_t& level ) const {
    detail::maximum_and_location( functionspace, field, maximum, glb_idx, level );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::minimumAndLocation( const Field& field, Vector& minimum,
                                                                      std::vector<gidx_t>& glb_idx ) const {
    detail::minimum_and_location( functionspace, field, minimum, glb_idx );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::maximumAndLocation( const Field& field, Vector& maximum,
                                                                      std::vector<gidx_t>& glb_idx ) const {
    detail::maximum_and_location( functionspace, field, maximum, glb_idx );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::minimumAndLocation( const Field& field, Vector& minimum,
                                                                      std::vector<gidx_t>& glb_idx,
                                                                      std::vector<idx_t>& level ) const {
    detail::minimum_and_location( functionspace, field, minimum, glb_idx, level );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::maximumAndLocation( const Field& field, Vector& maximum,
                                                                      std::vector<gidx_t>& glb_idx,
                                                                      std::vector<idx_t>& level ) const {
    detail::maximum_and_location( functionspace, field, maximum, glb_idx, level );
}

void NodeColumns::FieldStatistics::minimumAndLocationPerLevel( const Field& field, Field& column,
                                                               Field& glb_idx ) const {
    detail::minimum_and_location_per_level( functionspace, field, column, glb_idx );
}

void NodeColumns::FieldStatistics::maximumAndLocationPerLevel( const Field& field, Field& column,
                                                               Field& glb_idx ) const {
    detail::maximum_and_location_per_level( functionspace, field, column, glb_idx );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::mean( const Field& field, Value& mean, idx_t& N ) const {
    detail::mean( functionspace, field, mean, N );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::mean( const Field& field, Vector& mean, idx_t& N ) const {
    detail::mean( functionspace, field, mean, N );
}

void NodeColumns::FieldStatistics::meanPerLevel( const Field& field, Field& mean, idx_t& N ) const {
    detail::mean_per_level( functionspace, field, mean, N );
}

template <typename Value>
void NodeColumns::FieldStatisticsT<Value>::meanAndStandardDeviation( const Field& field, Value& mean, Value& stddev,
                                                                     idx_t& N ) const {
    detail::mean_and_standard_deviation( functionspace, field, mean, stddev, N );
}

template <typename Vector>
void NodeColumns::FieldStatisticsVectorT<Vector>::meanAndStandardDeviation( const Field& field, Vector& mean,
                                                                            Vector& stddev, idx_t& N ) const {
    detail::mean_and_standard_deviation( functionspace, field, mean, stddev, N );
}

void NodeColumns::FieldStatistics::meanAndStandardDeviationPerLevel( const Field& field, Field& mean, Field& stddev,
                                                                     idx_t& N ) const {
    detail::mean_and_standard_deviation_per_level( functionspace, field, mean, stddev, N );
}

template class NodeColumns::FieldStatisticsT<int>;
template class NodeColumns::FieldStatisticsT<long>;
template class NodeColumns::FieldStatisticsT<float>;
template class NodeColumns::FieldStatisticsT<double>;
// template class NodeColumns::FieldStatisticsT<unsigned long>;
template class NodeColumns::FieldStatisticsVectorT<std::vector<int>>;
template class NodeColumns::FieldStatisticsVectorT<std::vector<long>>;
template class NodeColumns::FieldStatisticsVectorT<std::vector<float>>;
template class NodeColumns::FieldStatisticsVectorT<std::vector<double>>;
// template class NodeColumns::FieldStatisticsVectorT< std::vector<unsigned
// long> >;

}  // namespace detail

NodeColumns::NodeColumns() : FunctionSpace(), functionspace_( nullptr ) {}

NodeColumns::NodeColumns( const FunctionSpace& functionspace ) :
    FunctionSpace( functionspace ),
    functionspace_( dynamic_cast<const detail::NodeColumns*>( get() ) ) {}

namespace {
detail::NodeColumns* make_functionspace( Mesh mesh, const eckit::Configuration& config ) {
    return new detail::NodeColumns( mesh, config );
}
}  // namespace

NodeColumns::NodeColumns( Mesh mesh ) :
    FunctionSpace( make_functionspace( mesh, util::NoConfig() ) ),
    functionspace_( dynamic_cast<const detail::NodeColumns*>( get() ) ) {}

NodeColumns::NodeColumns( Mesh mesh, const eckit::Configuration& config ) :
    FunctionSpace( make_functionspace( mesh, config ) ),
    functionspace_( dynamic_cast<const detail::NodeColumns*>( get() ) ) {}

idx_t NodeColumns::nb_nodes() const {
    return functionspace_->nb_nodes();
}

idx_t NodeColumns::nb_nodes_global() const {  // All MPI ranks will have same output
    return functionspace_->nb_nodes_global();
}

const Mesh& NodeColumns::mesh() const {
    return functionspace_->mesh();
}

mesh::Nodes& NodeColumns::nodes() const {
    return functionspace_->nodes();
}

// -- Parallelisation aware methods

const mesh::Halo& NodeColumns::halo() const {
    return functionspace_->halo();
}

void NodeColumns::haloExchange( FieldSet& fieldset, bool on_device ) const {
    functionspace_->haloExchange( fieldset, on_device );
}

void NodeColumns::haloExchange( Field& field, bool on_device ) const {
    functionspace_->haloExchange( field, on_device );
}

const parallel::HaloExchange& NodeColumns::halo_exchange() const {
    return functionspace_->halo_exchange();
}

void NodeColumns::gather( const FieldSet& local, FieldSet& global ) const {
    functionspace_->gather( local, global );
}

void NodeColumns::gather( const Field& local, Field& global ) const {
    functionspace_->gather( local, global );
}

const parallel::GatherScatter& NodeColumns::gather() const {
    return functionspace_->gather();
}

void NodeColumns::scatter( const FieldSet& global, FieldSet& local ) const {
    functionspace_->scatter( global, local );
}

void NodeColumns::scatter( const Field& global, Field& local ) const {
    functionspace_->scatter( global, local );
}

const parallel::GatherScatter& NodeColumns::scatter() const {
    return functionspace_->scatter();
}

std::string NodeColumns::checksum( const FieldSet& fieldset ) const {
    return functionspace_->checksum( fieldset );
}

std::string NodeColumns::checksum( const Field& field ) const {
    return functionspace_->checksum( field );
}

const parallel::Checksum& NodeColumns::checksum() const {
    return functionspace_->checksum();
}

}  // namespace functionspace
}  // namespace atlas
