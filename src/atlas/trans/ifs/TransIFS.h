/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "transi/trans.h"

#include "eckit/filesystem/PathName.h"

#include "atlas/array/LocalView.h"
#include "atlas/grid/Grid.h"
#include "atlas/trans/Trans.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field
}  // namespace atlas

namespace atlas {
namespace array {
class Array;
}
}  // namespace atlas

namespace atlas {
namespace functionspace {
class StructuredColumns;
class NodeColumns;
class Spectral;
}  // namespace functionspace
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {
class NodeColumns;
class Spectral;
}  // namespace detail
}  // namespace functionspace
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {
class TransPartitioner;
}
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace trans {
namespace detail {
void Assert( int code, const char* msg, const char* file, int line, const char* func );
#define ATLAS_TRANS_ASSERT( a ) atlas::trans::detail::Assert( !( a ), #a, __FILE__, __LINE__, __func__ )
}  // namespace detail
}  // namespace trans
}  // namespace atlas


//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class TransIFS : public trans::TransImpl {
private:
    typedef struct ::Trans_t Trans_t;

public:
    TransIFS( const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig() );
    TransIFS( const Grid&, const Domain&, const long truncation, const eckit::Configuration& = util::NoConfig() );
    TransIFS( const Cache&, const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig() );
    TransIFS( const Cache&, const Grid&, const Domain&, const long truncation,
              const eckit::Configuration& = util::NoConfig() );

    virtual ~TransIFS() override;
    operator ::Trans_t*() const { return trans(); }
    ::Trans_t* trans() const { return trans_.get(); }

    virtual int truncation() const override { return std::max( 0, trans_->nsmax ); }
    virtual size_t spectralCoefficients() const override { return trans_->nspec2g; }

    virtual const Grid& grid() const override { return grid_; }

    // pure virtual interface

    virtual void dirtrans( const Field& gpfield, Field& spfield,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans( const FieldSet& gpfields, FieldSet& spfields,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                       const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans( const Field& spfield, Field& gpfield,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans( const FieldSet& spfields, FieldSet& gpfields,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans_grad( const Field& spfield, Field& gradfield,
                                const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans_grad( const FieldSet& spfields, FieldSet& gradfields,
                                const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                       const eckit::Configuration& = util::NoConfig() ) const override;

    // -- IFS style API --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
 * @brief invtrans
 * @param nb_scalar_fields
 * @param scalar_spectra       [NSPEC2][nb_scalar_fields]
 * @param nb_vordiv_fields
 * @param vorticity_spectra    [NSPEC2][nb_vordiv_fields]
 * @param divergence_spectra   [NSPEC2][nb_vordiv_fields]
 * @param gp_fields  Ordering: [NGPBLKS][NFLD][NPROMA] if distributed,
 *                             [NFLD][NGPTOTG] if global ( add option::global()
 * )
 */
    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                           const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    /*!
 * @brief invtrans
 * @param nb_fields
 * @param scalar_spectra
 * @param scalar_fields
 */
    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    /*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
    virtual void invtrans( const int nb_vordiv_fields, const double vorticity_spectra[],
                           const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    /*!
 * @brief Direct transform of scalar fields
 */
    virtual void dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    /*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
    virtual void dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                           double divergence_spectra[], const eckit::Configuration& = util::NoConfig() ) const override;

    // implementations

public:
    void __dirtrans( const functionspace::StructuredColumns&, const Field& gpfield, const functionspace::Spectral&,
                     Field& spfield, const eckit::Configuration& = util::NoConfig() ) const;
    void __dirtrans( const functionspace::StructuredColumns&, const FieldSet& gpfields, const functionspace::Spectral&,
                     FieldSet& spfields, const eckit::Configuration& = util::NoConfig() ) const;

    void __dirtrans( const functionspace::NodeColumns&, const Field& gpfield, const functionspace::Spectral&,
                     Field& spfield, const eckit::Configuration& = util::NoConfig() ) const;
    void __dirtrans( const functionspace::NodeColumns&, const FieldSet& gpfields, const functionspace::Spectral&,
                     FieldSet& spfields, const eckit::Configuration& = util::NoConfig() ) const;

    void __dirtrans_wind2vordiv( const functionspace::NodeColumns&, const Field& gpwind, const functionspace::Spectral&,
                                 Field& spvor, Field& spdiv, const eckit::Configuration& = util::NoConfig() ) const;

    void __invtrans( const functionspace::Spectral&, const Field& spfield, const functionspace::StructuredColumns&,
                     Field& gpfield, const eckit::Configuration& = util::NoConfig() ) const;
    void __invtrans( const functionspace::Spectral&, const FieldSet& spfields, const functionspace::StructuredColumns&,
                     FieldSet& gpfields, const eckit::Configuration& = util::NoConfig() ) const;

    void __invtrans( const functionspace::Spectral&, const Field& spfield, const functionspace::NodeColumns&,
                     Field& gpfield, const eckit::Configuration& = util::NoConfig() ) const;
    void __invtrans( const functionspace::Spectral&, const FieldSet& spfields, const functionspace::NodeColumns&,
                     FieldSet& gpfields, const eckit::Configuration& = util::NoConfig() ) const;

    void __invtrans_vordiv2wind( const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                                 const functionspace::NodeColumns&, Field& gpwind,
                                 const eckit::Configuration& = util::NoConfig() ) const;

    void __invtrans_grad( const functionspace::Spectral& sp, const Field& spfield, const functionspace::NodeColumns& gp,
                          Field& gradfield, const eckit::Configuration& = util::NoConfig() ) const;

    void __invtrans_grad( const functionspace::Spectral& sp, const FieldSet& spfields,
                          const functionspace::NodeColumns& gp, FieldSet& gradfields,
                          const eckit::Configuration& = util::NoConfig() ) const;

public:
    void specnorm( const int nb_fields, const double spectra[], double norms[], int rank = 0 ) const;

protected:
    void assertCompatibleDistributions( const FunctionSpace& gp, const FunctionSpace& /*sp*/ ) const;

private:
    void ctor( const Grid&, long nsmax, const eckit::Configuration& );

    void ctor_rgg( const long nlat, const idx_t pl[], long nsmax, const eckit::Configuration& );

    void ctor_lonlat( const long nlon, const long nlat, long nsmax, const eckit::Configuration& );

    void ctor_spectral_only( long truncation, const eckit::Configuration& );

private:
    friend class grid::detail::partitioner::TransPartitioner;

    /// @brief Constructor for grid-only setup (e.g. for
    /// partitioning/parallelisation routines)
    TransIFS( const Grid& g, const eckit::Configuration& = util::NoConfig() );

    int handle() const { return trans_->handle; }
    int ndgl() const { return trans_->ndgl; }
    int nsmax() const { return trans_->nsmax; }
    int ngptot() const { return trans_->ngptot; }
    int ngptotg() const { return trans_->ngptotg; }
    int ngptotmx() const { return trans_->ngptotmx; }
    int nspec() const { return trans_->nspec; }
    int nspec2() const { return trans_->nspec2; }
    int nspec2g() const { return trans_->nspec2g; }
    int nspec2mx() const { return trans_->nspec2mx; }
    int n_regions_NS() const { return trans_->n_regions_NS; }
    int n_regions_EW() const { return trans_->n_regions_EW; }
    int nump() const { return trans_->nump; }
    int nproc() const { return trans_->nproc; }
    int myproc( int proc0 = 0 ) const { return trans_->myproc - 1 + proc0; }

    const int* nloen( int& size ) const {
        size = trans_->ndgl;
        ATLAS_TRANS_ASSERT( trans_->nloen != nullptr );
        return trans_->nloen;
    }

    array::LocalView<int, 1> nloen() const {
        ATLAS_TRANS_ASSERT( trans_->nloen != nullptr );
        return array::LocalView<int, 1>( trans_->nloen, array::make_shape( trans_->ndgl ) );
    }

    const int* n_regions( int& size ) const {
        size = trans_->n_regions_NS;
        ATLAS_TRANS_ASSERT( trans_->n_regions != nullptr );
        return trans_->n_regions;
    }

    array::LocalView<int, 1> n_regions() const {
        ATLAS_TRANS_ASSERT( trans_->n_regions != nullptr );
        return array::LocalView<int, 1>( trans_->n_regions, array::make_shape( trans_->n_regions_NS ) );
    }

    const int* nfrstlat( int& size ) const {
        size = trans_->n_regions_NS;
        if ( trans_->nfrstlat == nullptr ) ::trans_inquire( trans_.get(), "nfrstlat" );
        return trans_->nfrstlat;
    }

    array::LocalView<int, 1> nfrstlat() const {
        if ( trans_->nfrstlat == nullptr ) ::trans_inquire( trans_.get(), "nfrstlat" );
        return array::LocalView<int, 1>( trans_->nfrstlat, array::make_shape( trans_->n_regions_NS ) );
    }

    const int* nlstlat( int& size ) const {
        size = trans_->n_regions_NS;
        if ( trans_->nlstlat == nullptr ) ::trans_inquire( trans_.get(), "nlstlat" );
        return trans_->nlstlat;
    }

    array::LocalView<int, 1> nlstlat() const {
        if ( trans_->nlstlat == nullptr ) ::trans_inquire( trans_.get(), "nlstlat" );
        return array::LocalView<int, 1>( trans_->nlstlat, array::make_shape( trans_->n_regions_NS ) );
    }

    const int* nptrfrstlat( int& size ) const {
        size = trans_->n_regions_NS;
        if ( trans_->nptrfrstlat == nullptr ) ::trans_inquire( trans_.get(), "nptrfrstlat" );
        return trans_->nptrfrstlat;
    }

    array::LocalView<int, 1> nptrfrstlat() const {
        if ( trans_->nptrfrstlat == nullptr ) ::trans_inquire( trans_.get(), "nptrfrstlat" );
        return array::LocalView<int, 1>( trans_->nptrfrstlat, array::make_shape( trans_->n_regions_NS ) );
    }

    const int* nsta( int& sizef2, int& sizef1 ) const {
        sizef1 = trans_->ndgl + trans_->n_regions_NS - 1;
        sizef2 = trans_->n_regions_EW;
        if ( trans_->nsta == nullptr ) ::trans_inquire( trans_.get(), "nsta" );
        return trans_->nsta;
    }

    array::LocalView<int, 2> nsta() const {
        if ( trans_->nsta == nullptr ) ::trans_inquire( trans_.get(), "nsta" );
        return array::LocalView<int, 2>(
            trans_->nsta, array::make_shape( trans_->n_regions_EW, trans_->ndgl + trans_->n_regions_NS - 1 ) );
    }

    const int* nonl( int& sizef2, int& sizef1 ) const {
        sizef1 = trans_->ndgl + trans_->n_regions_NS - 1;
        sizef2 = trans_->n_regions_EW;
        if ( trans_->nonl == nullptr ) ::trans_inquire( trans_.get(), "nonl" );
        return trans_->nonl;
    }

    array::LocalView<int, 2> nonl() const {
        if ( trans_->nonl == nullptr ) ::trans_inquire( trans_.get(), "nonl" );
        return array::LocalView<int, 2>(
            trans_->nonl, array::make_shape( trans_->n_regions_EW, trans_->ndgl + trans_->n_regions_NS - 1 ) );
    }

    const int* nmyms( int& size ) const {
        size = trans_->nump;
        if ( trans_->nmyms == nullptr ) ::trans_inquire( trans_.get(), "nmyms" );
        return trans_->nmyms;
    }

    array::LocalView<int, 1> nmyms() const {
        if ( trans_->nmyms == nullptr ) ::trans_inquire( trans_.get(), "nmyms" );
        return array::LocalView<int, 1>( trans_->nmyms, array::make_shape( trans_->nump ) );
    }

    const int* nasm0( int& size ) const {
        size = trans_->nsmax + 1;  // +1 because zeroth wave included
        if ( trans_->nasm0 == nullptr ) ::trans_inquire( trans_.get(), "nasm0" );
        return trans_->nasm0;
    }

    array::LocalView<int, 1> nasm0() const {
        if ( trans_->nasm0 == nullptr ) ::trans_inquire( trans_.get(), "nasm0" );
        return array::LocalView<int, 1>( trans_->nasm0, array::make_shape( trans_->nsmax + 1 ) );
    }

    const int* nvalue( int& size ) const {
        size = trans_->nspec2;
        if ( trans_->nvalue == nullptr ) ::trans_inquire( trans_.get(), "nvalue" );
        return trans_->nvalue;
    }

    array::LocalView<int, 1> nvalue() const {
        if ( trans_->nvalue == nullptr ) ::trans_inquire( trans_.get(), "nvalue" );
        return array::LocalView<int, 1>( trans_->nvalue, array::make_shape( trans_->nspec2 ) );
    }

public:
    /*!
   * @brief distspec
   * @param nb_fields
   * @param origin
   * @param global_spectra
   * @param spectra
   */
    void distspec( const int nb_fields, const int origin[], const double global_spectra[], double spectra[] ) const;

    /*!
   * @brief gathspec
   * @param nb_fields
   * @param destination
   * @param spectra
   * @param global_spectra
   */
    void gathspec( const int nb_fields, const int destination[], const double spectra[],
                   double global_spectra[] ) const;

    /*!
   * @brief distgrid
   * @param nb_fields
   * @param origin
   * @param global_fields
   * @param fields
   */
    void distgrid( const int nb_fields, const int origin[], const double global_fields[], double fields[] ) const;

    /*!
   * @brief gathgrid
   * @param nb_fields
   * @param destination
   * @param fields
   * @param global_fields
   */
    void gathgrid( const int nb_fields, const int destination[], const double fields[], double global_fields[] ) const;

private:
    friend class functionspace::detail::Spectral;
    mutable std::shared_ptr<::Trans_t> trans_;
    grid::StructuredGrid grid_;
    const void* cache_{nullptr};
    size_t cachesize_{0};
};

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

extern "C" {
TransIFS* atlas__Trans__new( const Grid::Implementation* grid, int nsmax );
void atlas__Trans__delete( TransIFS* trans );
void atlas__Trans__distspec( const TransIFS* t, int nb_fields, int origin[], double global_spectra[],
                             double spectra[] );
void atlas__Trans__gathspec( const TransIFS* t, int nb_fields, int destination[], double spectra[],
                             double global_spectra[] );
void atlas__Trans__distgrid( const TransIFS* t, int nb_fields, int origin[], double global_fields[], double fields[] );
void atlas__Trans__gathgrid( const TransIFS* t, int nb_fields, int destination[], double fields[],
                             double global_fields[] );
void atlas__Trans__invtrans( const TransIFS* t, int nb_scalar_fields, double scalar_spectra[], int nb_vordiv_fields,
                             double vorticity_spectra[], double divergence_spectra[], double gp_fields[],
                             const eckit::Configuration* parameters );
void atlas__Trans__invtrans_scalar( const TransIFS* t, int nb_fields, double scalar_spectra[], double scalar_fields[] );
void atlas__Trans__invtrans_vordiv2wind( const TransIFS* t, int nb_fields, double vorticity_spectra[],
                                         double divergence_spectra[], double wind_fields[] );
void atlas__Trans__dirtrans_scalar( const TransIFS* t, int nb_fields, double scalar_fields[], double scalar_spectra[] );
void atlas__Trans__dirtrans_wind2vordiv( const TransIFS* t, int nb_fields, double wind_fields[],
                                         double vorticity_spectra[], double divergence_spectra[] );
void atlas__Trans__dirtrans_wind2vordiv_field( const TransIFS* This, const field::FieldImpl* gpwind,
                                               field::FieldImpl* spvor, field::FieldImpl* spdiv,
                                               const eckit::Configuration* parameters );
void atlas__Trans__specnorm( const TransIFS* t, int nb_fields, double spectra[], double norms[], int rank );
void atlas__Trans__dirtrans_fieldset( const TransIFS* This, const field::FieldSetImpl* gpfields,
                                      field::FieldSetImpl* spfields, const eckit::Configuration* parameters );
void atlas__Trans__dirtrans_field( const TransIFS* This, const field::FieldImpl* gpfield, field::FieldImpl* spfield,
                                   const eckit::Configuration* parameters );
void atlas__Trans__invtrans_fieldset( const TransIFS* This, const field::FieldSetImpl* spfields,
                                      field::FieldSetImpl* gpfields, const eckit::Configuration* parameters );
void atlas__Trans__invtrans_field( const TransIFS* This, const field::FieldImpl* spfield, field::FieldImpl* gpfield,
                                   const eckit::Configuration* parameters );
void atlas__Trans__invtrans_grad_field( const TransIFS* This, const field::FieldImpl* spfield,
                                        field::FieldImpl* gpfield, const eckit::Configuration* parameters );
void atlas__Trans__invtrans_vordiv2wind_field( const TransIFS* This, const field::FieldImpl* spvor,
                                               const field::FieldImpl* spdiv, field::FieldImpl* gpwind,
                                               const eckit::Configuration* parameters );

int atlas__Trans__handle( const TransIFS* trans );
int atlas__Trans__truncation( const TransIFS* This );
int atlas__Trans__nspec2( const TransIFS* This );
int atlas__Trans__nspec2g( const TransIFS* This );
int atlas__Trans__ngptot( const TransIFS* This );
int atlas__Trans__ngptotg( const TransIFS* This );
const Grid::Implementation* atlas__Trans__grid( const TransIFS* This );
}

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
