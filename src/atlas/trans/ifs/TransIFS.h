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

#include "eckit/filesystem/PathName.h"

#include "atlas/array/LocalView.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/trans/detail/TransImpl.h"

//-----------------------------------------------------------------------------
// Forward declarations

struct Trans_t;  // declared in "ectrans/transi.h"

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


//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class TransIFS : public trans::TransImpl {
private:
    typedef struct ::Trans_t Trans_t;

public:
    TransIFS(const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransIFS(const Grid&, const Domain&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransIFS(const Cache&, const Grid&, const long truncation, const eckit::Configuration& = util::NoConfig());
    TransIFS(const Cache&, const Grid&, const Domain&, const long truncation,
             const eckit::Configuration& = util::NoConfig());

    virtual ~TransIFS() override;
    operator ::Trans_t*() const { return trans(); }
    ::Trans_t* trans() const { return trans_.get(); }

    virtual int truncation() const override;
    virtual size_t nb_spectral_coefficients() const override;
    virtual size_t nb_spectral_coefficients_global() const override;

    virtual const Grid& grid() const override { return grid_; }
    virtual const functionspace::Spectral& spectral() const override;

    // pure virtual interface

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
    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                          const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    /*!
 * @brief invtrans
 * @param nb_fields
 * @param scalar_spectra
 * @param scalar_fields
 */
    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    /*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
    virtual void invtrans(const int nb_vordiv_fields, const double vorticity_spectra[],
                          const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const override;


    /*!
 * @brief invtrans_adj
 * @param nb_scalar_fields
 * @param scalar_spectra       [NSPEC2][nb_scalar_fields]
 * @param nb_vordiv_fields
 * @param vorticity_spectra    [NSPEC2][nb_vordiv_fields]
 * @param divergence_spectra   [NSPEC2][nb_vordiv_fields]
 * @param gp_fields  Ordering: [NGPBLKS][NFLD][NPROMA] if distributed,
 *                             [NFLD][NGPTOTG] if global ( add option::global()
 * )
 */
    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                              double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;


    /*!
 * @brief invtrans_adj
 * @param nb_fields
 * @param scalar_spectra
 * @param scalar_fields
 */
    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;


    /*!
 * @brief invtrans_adj
 * @param nb_scalar_fields
 * @param scalar_spectra       [NSPEC2][nb_scalar_fields]
 * @param nb_vordiv_fields
 * @param vorticity_spectra    [NSPEC2][nb_vordiv_fields]
 * @param divergence_spectra   [NSPEC2][nb_vordiv_fields]
 * @param gp_fields  Ordering: [NGPBLKS][NFLD][NPROMA] if distributed,
 *                             [NFLD][NGPTOTG] if global ( add option::global()
 * )
 */
    virtual void invtrans_adj(const int nb_vordiv_fields, const double gp_fields[], double vorticity_spectra[],
                              double divergence_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const override;

    /*!
 * @brief Direct transform of scalar fields
 */
    virtual void dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                          const eckit::Configuration& = util::NoConfig()) const override;

    /*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
    virtual void dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                          double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const override;

    // implementations

public:
    void __dirtrans(const functionspace::StructuredColumns&, const Field& gpfield, const functionspace::Spectral&,
                    Field& spfield, const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans(const functionspace::NodeColumns&, const Field& gpfield, const functionspace::Spectral&,
                    Field& spfield, const eckit::Configuration& = util::NoConfig()) const;

    void __dirtrans(const functionspace::StructuredColumns&, const FieldSet& gpfields, const functionspace::Spectral&,
                    FieldSet& spfields, const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans(const functionspace::NodeColumns&, const FieldSet& gpfields, const functionspace::Spectral&,
                    FieldSet& spfields, const eckit::Configuration& = util::NoConfig()) const;

    void __dirtrans_wind2vordiv(const functionspace::StructuredColumns&, const Field& gpwind,
                                const functionspace::Spectral&, Field& spvor, Field& spdiv,
                                const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans_wind2vordiv(const functionspace::NodeColumns&, const Field& gpwind, const functionspace::Spectral&,
                                Field& spvor, Field& spdiv, const eckit::Configuration& = util::NoConfig()) const;

    void __dirtrans_adj(const functionspace::Spectral&, const Field& spfield,
                        const functionspace::StructuredColumns&, Field& gpfield,
                        const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans_adj(const functionspace::Spectral&, const Field& spfield,
                        const functionspace::NodeColumns&, Field& gpfield,
                        const eckit::Configuration& = util::NoConfig()) const;

    void __dirtrans_adj(const functionspace::Spectral&, const FieldSet& spfields,
                        const functionspace::StructuredColumns&, FieldSet& gpfields,
                        const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans_adj(const functionspace::Spectral&, const FieldSet& spfields,
                        const functionspace::NodeColumns&, FieldSet& gpfields,
                        const eckit::Configuration& = util::NoConfig()) const;

    void __dirtrans_wind2vordiv_adj(const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                                    const functionspace::StructuredColumns&, Field& gpwind,
                                    const eckit::Configuration& = util::NoConfig()) const;
    void __dirtrans_wind2vordiv_adj(const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                                    const functionspace::NodeColumns&, Field& gpwind,
                                    const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans(const functionspace::Spectral&, const Field& spfield, const functionspace::StructuredColumns&,
                    Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans(const functionspace::Spectral&, const Field& spfield, const functionspace::NodeColumns&,
                    Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans(const functionspace::Spectral&, const FieldSet& spfields, const functionspace::StructuredColumns&,
                    FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans(const functionspace::Spectral&, const FieldSet& spfields, const functionspace::NodeColumns&,
                    FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_grad(const functionspace::Spectral& sp, const Field& spfield,
                         const functionspace::StructuredColumns& gp, Field& gradfield,
                         const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_grad(const functionspace::Spectral& sp, const Field& spfield, const functionspace::NodeColumns& gp,
                         Field& gradfield, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_grad(const functionspace::Spectral& sp, const FieldSet& spfields,
                         const functionspace::StructuredColumns& gp, FieldSet& gradfields,
                         const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_grad(const functionspace::Spectral& sp, const FieldSet& spfields,
                         const functionspace::NodeColumns& gp, FieldSet& gradfields,
                         const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_vordiv2wind(const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                                const functionspace::StructuredColumns&, Field& gpwind,
                                const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_vordiv2wind(const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                                const functionspace::NodeColumns&, Field& gpwind,
                                const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_adj(const functionspace::Spectral&, Field& spfield, const functionspace::StructuredColumns&,
                        const Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_adj(const functionspace::Spectral&, Field& spfield, const functionspace::NodeColumns&,
                        const Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_adj(const functionspace::Spectral&, FieldSet& spfields, const functionspace::StructuredColumns&,
                        const FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_adj(const functionspace::Spectral&, FieldSet& spfields, const functionspace::NodeColumns&,
                        const FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_grad_adj(const functionspace::Spectral& sp, Field& spfield,
                             const functionspace::StructuredColumns& gp, const Field& gradfield,
                             const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_grad_adj(const functionspace::Spectral& sp, Field& spfield, const functionspace::NodeColumns& gp,
                             const Field& gradfield, const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_grad_adj(const functionspace::Spectral& sp, FieldSet& spfields,
                             const functionspace::StructuredColumns& gp, const FieldSet& gradfields,
                             const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_grad_adj(const functionspace::Spectral& sp, FieldSet& spfields,
                             const functionspace::NodeColumns& gp, const FieldSet& gradfields,
                             const eckit::Configuration& = util::NoConfig()) const;

    void __invtrans_vordiv2wind_adj(const functionspace::Spectral&, Field& spvor, Field& spdiv,
                                    const functionspace::StructuredColumns&, const Field& gpwind,
                                    const eckit::Configuration& = util::NoConfig()) const;
    void __invtrans_vordiv2wind_adj(const functionspace::Spectral&, Field& spvor, Field& spdiv,
                                    const functionspace::NodeColumns&, const Field& gpwind,
                                    const eckit::Configuration& = util::NoConfig()) const;

public:
    void specnorm(const int nb_fields, const double spectra[], double norms[], int rank = 0) const;

protected:
    void assertCompatibleDistributions(const FunctionSpace& gp, const FunctionSpace& /*sp*/) const;

private:
    void ctor(const Grid&, long nsmax, const eckit::Configuration&);

    void ctor_rgg(const long nlat, const idx_t pl[], long nsmax, const eckit::Configuration&);

    void ctor_lonlat(const long nlon, const long nlat, long nsmax, const eckit::Configuration&);

private:
    friend class grid::detail::partitioner::TransPartitioner;

    /// @brief Constructor for grid-only setup (e.g. for
    /// partitioning/parallelisation routines)
    TransIFS(const Grid& g, const eckit::Configuration& = util::NoConfig());

    virtual int handle() const override;
    int ndgl() const;
    int nsmax() const;
    int ngptot() const;
    int ngptotg() const;
    int ngptotmx() const;
    int nspec() const;
    int nspec2() const;
    int nspec2g() const;
    int nspec2mx() const;
    int n_regions_NS() const;
    int n_regions_EW() const;
    int nump() const;
    int nproc() const;
    int myproc(int proc0 = 0) const;

    const int* nloen(int& size) const;

    array::LocalView<int, 1> nloen() const;

    const int* n_regions(int& size) const;

    array::LocalView<int, 1> n_regions() const;

    const int* nfrstlat(int& size) const;

    array::LocalView<int, 1> nfrstlat() const;

    const int* nlstlat(int& size) const;

    array::LocalView<int, 1> nlstlat() const;

    const int* nptrfrstlat(int& size) const;

    array::LocalView<int, 1> nptrfrstlat() const;

    const int* nsta(int& sizef2, int& sizef1) const;

    array::LocalView<int, 2> nsta() const;

    const int* nonl(int& sizef2, int& sizef1) const;

    array::LocalView<int, 2> nonl() const;

    const int* nmyms(int& size) const;

    array::LocalView<int, 1> nmyms() const;

    const int* nasm0(int& size) const;

    array::LocalView<int, 1> nasm0() const;

    const int* nvalue(int& size) const;

    array::LocalView<int, 1> nvalue() const;

public:
    /*!
   * @brief distspec
   * @param nb_fields
   * @param origin
   * @param global_spectra
   * @param spectra
   */
    void distspec(const int nb_fields, const int origin[], const double global_spectra[], double spectra[]) const;

    /*!
   * @brief gathspec
   * @param nb_fields
   * @param destination
   * @param spectra
   * @param global_spectra
   */
    void gathspec(const int nb_fields, const int destination[], const double spectra[], double global_spectra[]) const;

    /*!
   * @brief distgrid
   * @param nb_fields
   * @param origin
   * @param global_fields
   * @param fields
   */
    void distgrid(const int nb_fields, const int origin[], const double global_fields[], double fields[]) const;

    /*!
   * @brief gathgrid
   * @param nb_fields
   * @param destination
   * @param fields
   * @param global_fields
   */
    void gathgrid(const int nb_fields, const int destination[], const double fields[], double global_fields[]) const;

private:
    friend class functionspace::detail::Spectral;
    mutable std::shared_ptr<::Trans_t> trans_;
    StructuredGrid grid_;
    const void* cache_{nullptr};
    size_t cachesize_{0};

protected:
    mutable functionspace::Spectral spectral_;
};

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
