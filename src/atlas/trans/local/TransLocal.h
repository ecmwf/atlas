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

#include <memory>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/grid/Grid.h"
#include "atlas/linalg/dense/Backend.h"
#include "atlas/trans/detail/TransImpl.h"

#define TRANSLOCAL_DGEMM2 0

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
class StructuredGrid;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

namespace detail {
struct FFTW_Data;
}

class LegendreCacheCreatorLocal;
int fourier_truncation(const int truncation,  // truncation
                       const int nx,          // number of longitudes
                       const int nxmax,       // maximum nx
                       const int ndgl,        // number of latitudes
                       const double lat,      // latitude in radian
                       const bool regular);   // regular grid

//-----------------------------------------------------------------------------

/// @class TransLocal
///
/// Local spherical harmonics transformations to any grid
/// Optimisations are present for structured grids
/// For global grids, please consider using TransIFS instead.
///
/// @todo:
///  - support multiple fields
///  - support atlas::Field and atlas::FieldSet based on function spaces
///
/// @note: Direct transforms are not implemented and cannot be unless
///        the grid is global. There are no plans to support this at the moment.
///
/// @note: The matrix_multiply (GEMM) implementation can be configured within the Configuration argument in the constructor
///        using "matrix_multiply" key or if not given, it will use the atlas::linalg::dense::current_backend(),
///        evaluated at invocation time. To reset the current_backend at any time:
///             atlas::linalg::dense::current_backend( type );
///        Possible values for type or "matrix_multiply" are e.g.
///        - "eckit_linalg" which uses the configured default backend of eckit::linalg::LinearAlgebra
///        - "generic" : "generic" backend for eckit::linalg::LinearAlgebra
///        - "lapack"  : "lapack"  backend for eckit::linalg::LinearAlgebra
///        - "openmp"  : "openmp"  backend for eckit::linalg::LinearAlgebra, or "generic" if "openmp" is not available.
///        - "eigen"   : "eigen"   backend for eckit::linalg::LinearAlgebra

class TransLocal : public trans::TransImpl {
public:
    TransLocal(const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransLocal(const Grid&, const Domain&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransLocal(const Cache&, const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransLocal(const Cache&, const Grid&, const Domain&, const long truncation,
               const eckit::Configuration& = util::NoConfig());

    virtual ~TransLocal() override;

    virtual int truncation() const override { return truncation_; }

    virtual size_t nb_spectral_coefficients() const override { return (truncation_ + 1) * (truncation_ + 2); }
    virtual size_t nb_spectral_coefficients_global() const override { return (truncation_ + 1) * (truncation_ + 2); }

    virtual const Grid& grid() const override { return grid_; }
    virtual const functionspace::Spectral& spectral() const override;

    virtual void invtrans(const Field& spfield, Field& gpfield,
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans(const FieldSet& spfields, FieldSet& gpfields,
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_grad(const Field& spfield, Field& gradfield,
                               const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_grad(const FieldSet& spfields, FieldSet& gradfields,
                               const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_vordiv2wind(const Field& spvor, const Field& spdiv, Field& gpwind,
                                      const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_adj(const Field& gpfield, Field& spfield,
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_adj(const FieldSet& gpfields, FieldSet& spfields,
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_grad_adj(const Field& gradfield, Field& spfield,
                                   const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_grad_adj(const FieldSet& spfields, FieldSet& gradfields,
                                   const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_vordiv2wind_adj(const Field& gpwind, Field& spvor, Field& spdiv,
                                          const eckit::Configuration& = util::NoConfig()) const override;

    // -- IFS style API --

    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                          const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans(const int nb_vordiv_fields, const double vorticity_spectra[],
                          const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                              double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void invtrans_adj(const int nb_vordiv_fields, const double gp_fields[], double vorticity_spectra[],
                              double divergence_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;

    // -- NOT SUPPORTED -- //

    virtual void dirtrans(const Field& gpfield, Field& spfield,
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans(const FieldSet& gpfields, FieldSet& spfields,
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans_wind2vordiv(const Field& gpwind, Field& spvor, Field& spdiv,
                                      const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans_adj(const Field& spfield, Field& gpfield,
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans_adj(const FieldSet& spfields, FieldSet& gpfields,
                              const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans_wind2vordiv_adj(const Field& spvor, const Field& spdiv, Field& gpwind,
                                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    virtual void dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                          double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const override;

private:
    int posMethod(const int jfld, const int imag, const int jlat, const int jm, const int nb_fields,
                  const int nlats) const {
#if !TRANSLOCAL_DGEMM2
        return imag + 2 * (jm + (truncation_ + 1) * (jlat + nlats * jfld));
#else
        return jfld + nb_fields * (jlat + nlats * (imag + 2 * (jm)));
#endif
    }

    void invtrans_legendre(const int truncation, const int nlats, const int nb_fields, const int nb_vordiv_fields,
                           const double scalar_spectra[], double scl_fourier[],
                           const eckit::Configuration& config) const;

    void invtrans_fourier_regular(const int nlats, const int nlons, const int nb_fields, double scl_fourier[],
                                  double gp_fields[], const eckit::Configuration& config) const;

    void invtrans_fourier_reduced(const int nlats, const StructuredGrid& g, const int nb_fields, double scl_fourier[],
                                  double gp_fields[], const eckit::Configuration& config) const;

    void invtrans_unstructured_precomp(const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                                       const double scalar_spectra[], double gp_fields[],
                                       const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_unstructured(const int truncation, const int nb_fields, const int nb_vordiv_fields,
                               const double scalar_spectra[], double gp_fields[],
                               const eckit::Configuration& config) const;

    void invtrans_uv(const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                     const double scalar_spectra[], double gp_fields[],
                     const eckit::Configuration& = util::NoConfig()) const;

    bool warning(const eckit::Configuration& = util::NoConfig()) const;

    friend class LegendreCacheCreatorLocal;

private:
    mutable functionspace::Spectral spectral_;
    Grid grid_;
    Grid gridGlobal_;
    bool useFFT_;
    bool dgemmMethod1_;
    bool unstruct_precomp_;
    bool no_symmetry_;
    int truncation_;
    idx_t nlatsNH_;
    idx_t nlatsSH_;
    idx_t nlatsLeg_;
    idx_t nlatsLegReduced_;
    idx_t nlatsLegDomain_;
    std::vector<idx_t> jlonMin_;
    idx_t jlatMin_;
    idx_t jlatMinLeg_;
    int nlonsMaxGlobal_;
    std::vector<idx_t> nlonsGlobal_;
    std::vector<idx_t> nlat0_;
    idx_t nlatsGlobal_;
    bool precompute_;
    double* legendre_;
    double* legendre_sym_;
    double* legendre_asym_;
    double* fourier_;
    double* fouriertp_;
    std::vector<size_t> legendre_begin_;
    std::vector<size_t> legendre_sym_begin_;
    std::vector<size_t> legendre_asym_begin_;

    Cache cache_;
    Cache export_legendre_;
    const void* legendre_cache_{nullptr};
    size_t legendre_cachesize_{0};
    const void* fft_cache_{nullptr};
    size_t fft_cachesize_{0};

    std::unique_ptr<detail::FFTW_Data> fftw_;

    std::string linalg_backend_;
    int warning_ = 0;
};

//-----------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
