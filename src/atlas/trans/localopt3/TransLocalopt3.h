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

#include <vector>

#include "atlas/array.h"
#include "atlas/grid/Grid.h"
#include "atlas/trans/Trans.h"
#if ATLAS_HAVE_FFTW
#include <fftw3.h>
#endif

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

/// @class TransLocalopt3
///
/// Localopt3 spherical harmonics transformations to any grid
/// Optimisations are present for structured grids
/// For global grids, please consider using TransIFS instead.
///
/// @todo:
///  - support multiple fields
///  - support atlas::Field and atlas::FieldSet based on function spaces
///
/// @note: Direct transforms are not implemented and cannot be unless
///        the grid is global. There are no plans to support this at the moment.
class TransLocalopt3 : public trans::TransImpl {
public:
    TransLocalopt3( const Grid& g, const long truncation, const eckit::Configuration& = util::NoConfig() );
    TransLocalopt3( const Cache&, const Grid& g, const long truncation,
                    const eckit::Configuration& = util::NoConfig() );

    virtual ~TransLocalopt3();

    virtual int truncation() const override { return truncation_; }
    virtual size_t spectralCoefficients() const override { return ( truncation_ + 1 ) * ( truncation_ + 2 ); }

    virtual const Grid& grid() const override { return grid_; }

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

    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                           const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void invtrans( const int nb_vordiv_fields, const double vorticity_spectra[],
                           const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    // -- NOT SUPPORTED -- //

    virtual void dirtrans( const Field& gpfield, Field& spfield,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans( const FieldSet& gpfields, FieldSet& spfields,
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                       const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                           const eckit::Configuration& = util::NoConfig() ) const override;

    virtual void dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                           double divergence_spectra[], const eckit::Configuration& = util::NoConfig() ) const override;

private:
    void invtrans_uv( const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                      const double scalar_spectra[], double gp_fields[],
                      const eckit::Configuration& = util::NoConfig() ) const;

private:
    Grid grid_;
    bool useFFT_;
    bool dgemmMethod1_;
    int truncation_;
    int nlatsNH_;
    int nlatsSH_;
    int nlatsLeg_;
    bool precompute_;
    double* legendre_sym_;
    double* legendre_asym_;
    double* fourier_;
    double* fouriertp_;
    std::vector<size_t> legendre_begin_;
    std::vector<size_t> legendre_sym_begin_;
    std::vector<size_t> legendre_asym_begin_;
#if ATLAS_HAVE_FFTW
    fftw_complex* fft_in_;
    double* fft_out_;
    fftw_plan plan_;
#endif
};

//-----------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
