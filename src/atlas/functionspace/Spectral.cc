/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/os/BackTrace.h"
#include "eckit/utils/MD5.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"


#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/TransIFS.h"
namespace {
void trans_check( const int code, const char* msg, const eckit::CodeLocation& location ) {
    if ( code != TRANS_SUCCESS ) {
        std::stringstream errmsg;
        errmsg << "atlas::trans ERROR: " << msg << " failed: \n";
        errmsg << ::trans_error_msg( code );
        atlas::throw_Exception( errmsg.str(), location );
    }
}
#define TRANS_CHECK( CALL ) trans_check( CALL, #CALL, Here() )
}  // namespace
#endif

namespace atlas {
namespace functionspace {
namespace detail {

#if ATLAS_HAVE_TRANS
class Spectral::Parallelisation {
public:
    Parallelisation( const std::shared_ptr<::Trans_t> other ) : trans_( other ) {}

    Parallelisation( int truncation ) {
        trans_ = std::shared_ptr<::Trans_t>( new ::Trans_t, [](::Trans_t* p ) {
            TRANS_CHECK(::trans_delete( p ) );
            delete p;
        } );
        TRANS_CHECK(::trans_new( trans_.get() ) );
        TRANS_CHECK(::trans_set_trunc( trans_.get(), truncation ) );
        TRANS_CHECK(::trans_use_mpi( mpi::comm().size() > 1 ) );
        TRANS_CHECK(::trans_setup( trans_.get() ) );
    }

    int nb_spectral_coefficients_global() const { return trans_->nspec2g; }
    int nb_spectral_coefficients() const { return trans_->nspec2; }

    std::string distribution() const { return "trans"; }
    operator ::Trans_t*() const { return trans_.get(); }
    std::shared_ptr<::Trans_t> trans_;
};
#else
class Spectral::Parallelisation {
public:
    Parallelisation( int truncation ) : truncation_( truncation ) {}
    int nb_spectral_coefficients_global() const { return ( truncation_ + 1 ) * ( truncation_ + 2 ); }
    int nb_spectral_coefficients() const { return nb_spectral_coefficients_global(); }
    int truncation_;
    std::string distribution() const { return "serial"; }
};
#endif

void Spectral::set_field_metadata( const eckit::Configuration& config, Field& field ) const {
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

    field.set_levels( config_levels( config ) );
    field.set_variables( 0 );
}

idx_t Spectral::config_size( const eckit::Configuration& config ) const {
    idx_t size = nb_spectral_coefficients();
    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            size = ( idx_t( mpi::comm().rank() ) == owner ? nb_spectral_coefficients_global() : 0 );
        }
    }
    return size;
}

// ----------------------------------------------------------------------

Spectral::Spectral( const eckit::Configuration& config ) :
    Spectral::Spectral( config.getUnsigned( "truncation" ), config ) {}

// ----------------------------------------------------------------------

Spectral::Spectral( const int truncation, const eckit::Configuration& config ) :
    nb_levels_( 0 ),
    truncation_( truncation ),
    parallelisation_( new Parallelisation( truncation_ ) ) {
    config.get( "levels", nb_levels_ );
}

Spectral::Spectral( const trans::Trans& trans, const eckit::Configuration& config ) :
    nb_levels_( 0 ),
    truncation_( trans.truncation() ),
#if ATLAS_HAVE_TRANS
    parallelisation_( new Parallelisation( dynamic_cast<const trans::TransIFS&>( *trans.get() ).trans_ ) ) {
#else
    parallelisation_( new Parallelisation( truncation_ ) ) {
#endif
    config.get( "levels", nb_levels_ );
}

Spectral::~Spectral() {}

std::string Spectral::distribution() const {
    return parallelisation_->distribution();
}

size_t Spectral::footprint() const {
    size_t size = sizeof( *this );
    // TODO
    return size;
}

idx_t Spectral::nb_spectral_coefficients() const {
    return parallelisation_->nb_spectral_coefficients();
}

idx_t Spectral::nb_spectral_coefficients_global() const {
    return parallelisation_->nb_spectral_coefficients_global();
}

array::DataType Spectral::config_datatype( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) throw_Exception( "datatype missing", Here() );
    return array::DataType( kind );
}

std::string Spectral::config_name( const eckit::Configuration& config ) const {
    std::string name;
    config.get( "name", name );
    return name;
}

idx_t Spectral::config_levels( const eckit::Configuration& config ) const {
    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    return levels;
}

Field Spectral::createField( const eckit::Configuration& options ) const {
    array::ArrayShape array_shape;

    idx_t nb_spec_coeffs = config_size( options );
    array_shape.push_back( nb_spec_coeffs );

    idx_t levels = config_levels( options );
    if ( levels ) array_shape.push_back( levels );

    Field field = Field( config_name( options ), config_datatype( options ), array_shape );

    set_field_metadata( options, field );
    return field;
}

Field Spectral::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | option::levels( other.levels() ) | config );
}

void Spectral::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const {
    ATLAS_ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& loc = local_fieldset[f];
        if ( loc.datatype() != array::DataType::str<double>() ) {
            std::stringstream err;
            err << "Cannot gather spectral field " << loc.name() << " of datatype " << loc.datatype().str() << ".";
            err << "Only " << array::DataType::str<double>() << " supported.";
            throw_Exception( err.str(), Here() );
        }

#if ATLAS_HAVE_TRANS
        Field& glb = global_fieldset[f];
        idx_t root = 0;
        idx_t rank = static_cast<idx_t>( mpi::comm().rank() );
        glb.metadata().get( "owner", root );
        ATLAS_ASSERT( loc.shape( 0 ) == nb_spectral_coefficients() );
        if ( rank == root ) ATLAS_ASSERT( glb.shape( 0 ) == nb_spectral_coefficients_global() );
        std::vector<int> nto( 1, root + 1 );
        if ( loc.rank() > 1 ) {
            nto.resize( loc.stride( 0 ) );
            for ( size_t i = 0; i < nto.size(); ++i )
                nto[i] = root + 1;
        }

        struct ::GathSpec_t args = new_gathspec( *parallelisation_ );
        args.nfld                = nto.size();
        args.rspecg              = glb.data<double>();
        args.nto                 = nto.data();
        args.rspec               = loc.data<double>();
        TRANS_CHECK(::trans_gathspec( &args ) );
#else

        throw_Exception(
            "Cannot gather spectral fields because Atlas has "
            "not been compiled with TRANS support." );
#endif
    }
}
void Spectral::gather( const Field& local, Field& global ) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add( local );
    global_fields.add( global );
    gather( local_fields, global_fields );
}

void Spectral::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const {
    ATLAS_ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& glb = global_fieldset[f];
        Field& loc       = local_fieldset[f];
        if ( loc.datatype() != array::DataType::str<double>() ) {
            std::stringstream err;
            err << "Cannot scatter spectral field " << glb.name() << " of datatype " << glb.datatype().str() << ".";
            err << "Only " << array::DataType::str<double>() << " supported.";
            throw_Exception( err.str(), Here() );
        }

#if ATLAS_HAVE_TRANS
        idx_t root = 0;
        idx_t rank = static_cast<idx_t>( mpi::comm().rank() );

        glb.metadata().get( "owner", root );
        ATLAS_ASSERT( loc.shape( 0 ) == nb_spectral_coefficients() );
        if ( rank == root ) ATLAS_ASSERT( glb.shape( 0 ) == nb_spectral_coefficients_global() );
        std::vector<int> nfrom( 1, root + 1 );
        if ( loc.rank() > 1 ) {
            nfrom.resize( loc.stride( 0 ) );
            for ( size_t i = 0; i < nfrom.size(); ++i )
                nfrom[i] = root + 1;
        }

        struct ::DistSpec_t args = new_distspec( *parallelisation_ );
        args.nfld                = int( nfrom.size() );
        args.rspecg              = glb.data<double>();
        args.nfrom               = nfrom.data();
        args.rspec               = loc.data<double>();
        TRANS_CHECK(::trans_distspec( &args ) );

        glb.metadata().broadcast( loc.metadata(), root );
        loc.metadata().set( "global", false );
#else
        throw_Exception(
            "Cannot scatter spectral fields because Atlas has "
            "not been compiled with TRANS support." );
#endif
    }
}
void Spectral::scatter( const Field& global, Field& local ) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add( global );
    local_fields.add( local );
    scatter( global_fields, local_fields );
}

std::string Spectral::checksum( const FieldSet& ) const {
    eckit::MD5 md5;
    ATLAS_NOTIMPLEMENTED;
}
std::string Spectral::checksum( const Field& field ) const {
    FieldSet fieldset;
    fieldset.add( field );
    return checksum( fieldset );
}

void Spectral::norm( const Field& field, double& norm, int rank ) const {
#if ATLAS_HAVE_TRANS
    ATLAS_ASSERT( std::min<int>( 1, field.levels() ) == 1 );
    struct ::SpecNorm_t args = new_specnorm( *parallelisation_ );
    args.nfld                = 1;
    args.rspec               = field.data<double>();
    args.rnorm               = &norm;
    args.nmaster             = rank + 1;
    TRANS_CHECK(::trans_specnorm( &args ) );
#else
    throw_Exception(
        "Cannot compute spectral norms because Atlas has not "
        "been compiled with TRANS support." );
#endif
}
void Spectral::norm( const Field& field, double norm_per_level[], int rank ) const {
#if ATLAS_HAVE_TRANS
    ATLAS_ASSERT( std::min<int>( 1, field.levels() ) == 1 );
    struct ::SpecNorm_t args = new_specnorm( *parallelisation_ );
    args.nfld                = std::max<int>( 1, field.levels() );
    args.rspec               = field.data<double>();
    args.rnorm               = norm_per_level;
    args.nmaster             = rank + 1;
    TRANS_CHECK(::trans_specnorm( &args ) );
#else
    throw_Exception(
        "Cannot compute spectral norms because Atlas has not "
        "been compiled with TRANS support." );
#endif
}
void Spectral::norm( const Field& field, std::vector<double>& norm_per_level, int rank ) const {
#if ATLAS_HAVE_TRANS
    norm_per_level.resize( std::min<int>( 1, field.levels() ) );
    struct ::SpecNorm_t args = new_specnorm( *parallelisation_ );
    args.nfld                = int( norm_per_level.size() );
    args.rspec               = field.data<double>();
    args.rnorm               = norm_per_level.data();
    args.nmaster             = rank + 1;
    TRANS_CHECK(::trans_specnorm( &args ) );
#else
    throw_Exception(
        "Cannot compute spectral norms because Atlas has not "
        "been compiled with TRANS support." );
#endif
}

}  // namespace detail

// ----------------------------------------------------------------------

Spectral::Spectral( const FunctionSpace& functionspace ) :
    FunctionSpace( functionspace ),
    functionspace_( dynamic_cast<const detail::Spectral*>( get() ) ) {}

Spectral::Spectral( const eckit::Configuration& config ) :
    FunctionSpace( new detail::Spectral( config ) ),
    functionspace_( dynamic_cast<const detail::Spectral*>( get() ) ) {}

Spectral::Spectral( const int truncation, const eckit::Configuration& config ) :
    FunctionSpace( new detail::Spectral( truncation, config ) ),
    functionspace_( dynamic_cast<const detail::Spectral*>( get() ) ) {}

Spectral::Spectral( const trans::Trans& trans, const eckit::Configuration& config ) :
    FunctionSpace( new detail::Spectral( trans, config ) ),
    functionspace_( dynamic_cast<const detail::Spectral*>( get() ) ) {}

idx_t Spectral::nb_spectral_coefficients() const {
    return functionspace_->nb_spectral_coefficients();
}

idx_t Spectral::nb_spectral_coefficients_global() const {
    return functionspace_->nb_spectral_coefficients_global();
}

int Spectral::truncation() const {
    return functionspace_->truncation();
}

void Spectral::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const {
    functionspace_->gather( local_fieldset, global_fieldset );
}

void Spectral::gather( const Field& local, Field& global ) const {
    functionspace_->gather( local, global );
}

void Spectral::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const {
    functionspace_->scatter( global_fieldset, local_fieldset );
}

void Spectral::scatter( const Field& global, Field& local ) const {
    functionspace_->scatter( global, local );
}

std::string Spectral::checksum( const FieldSet& fieldset ) const {
    return functionspace_->checksum( fieldset );
}

std::string Spectral::checksum( const Field& field ) const {
    return functionspace_->checksum( field );
}

void Spectral::norm( const Field& field, double& norm, int rank ) const {
    functionspace_->norm( field, norm, rank );
}

void Spectral::norm( const Field& field, double norm_per_level[], int rank ) const {
    functionspace_->norm( field, norm_per_level, rank );
}

void Spectral::norm( const Field& field, std::vector<double>& norm_per_level, int rank ) const {
    functionspace_->norm( field, norm_per_level, rank );
}

// ----------------------------------------------------------------------

extern "C" {
const detail::Spectral* atlas__SpectralFunctionSpace__new__config( const eckit::Configuration* config ) {
    ATLAS_ASSERT( config != nullptr );
    return new detail::Spectral( *config );
}

const detail::Spectral* atlas__SpectralFunctionSpace__new__trans( trans::TransImpl* trans,
                                                                  const eckit::Configuration* config ) {
    ATLAS_ASSERT( trans != nullptr );
    ATLAS_ASSERT( config != nullptr );
    return new detail::Spectral( trans::Trans( trans ), *config );
}

void atlas__SpectralFunctionSpace__delete( detail::Spectral* This ) {
    ATLAS_ASSERT( This != nullptr );
    delete This;
}

field::FieldImpl* atlas__fs__Spectral__create_field( const detail::Spectral* This,
                                                     const eckit::Configuration* options ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( options );
    field::FieldImpl* field;
    {
        Field f = This->createField( *options );
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

void atlas__SpectralFunctionSpace__gather( const detail::Spectral* This, const field::FieldImpl* local,
                                           field::FieldImpl* global ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( global != nullptr );
    ATLAS_ASSERT( local != nullptr );
    const Field l( local );
    Field g( global );
    This->gather( l, g );
}

void atlas__SpectralFunctionSpace__scatter( const detail::Spectral* This, const field::FieldImpl* global,
                                            field::FieldImpl* local ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( global != nullptr );
    ATLAS_ASSERT( local != nullptr );
    const Field g( global );
    Field l( local );
    This->scatter( g, l );
}

void atlas__SpectralFunctionSpace__gather_fieldset( const detail::Spectral* This, const field::FieldSetImpl* local,
                                                    field::FieldSetImpl* global ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( global != nullptr );
    ATLAS_ASSERT( local != nullptr );
    const FieldSet l( local );
    FieldSet g( global );
    This->gather( l, g );
}

void atlas__SpectralFunctionSpace__scatter_fieldset( const detail::Spectral* This, const field::FieldSetImpl* global,
                                                     field::FieldSetImpl* local ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( global != nullptr );
    ATLAS_ASSERT( local != nullptr );
    const FieldSet g( global );
    FieldSet l( local );
    This->scatter( g, l );
}

void atlas__SpectralFunctionSpace__norm( const detail::Spectral* This, const field::FieldImpl* field, double norm[],
                                         int rank ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( field != nullptr );
    ATLAS_ASSERT( norm != nullptr );
    This->norm( field, norm, rank );
}
}

}  // namespace functionspace
}  // namespace atlas
