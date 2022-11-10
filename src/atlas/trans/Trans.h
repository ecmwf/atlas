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

#include <iosfwd>

#include "atlas/library/config.h"
#include "atlas/trans/Cache.h"
#include "atlas/util/Config.h"
#include "atlas/util/ObjectHandle.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
class Domain;
namespace functionspace {
class Spectral;
}
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//----------------------------------------------------------------------------------------------------------------------

class TransImpl;
class Trans : DOXYGEN_HIDE(public util::ObjectHandle<TransImpl>) {
public:
    static void listBackends(std::ostream&);
    static bool hasBackend(const std::string&);
    static void backend(const std::string&);
    static std::string backend();
    static void config(const eckit::Configuration&);
    static const eckit::Configuration& config();

    using Handle::Handle;
    Trans() = default;

    Trans(const FunctionSpace& gp, const FunctionSpace& sp, const eckit::Configuration& = util::NoConfig());
    Trans(const Grid&, int truncation, const eckit::Configuration& = util::NoConfig());
    Trans(const Grid&, const Domain&, int truncation, const eckit::Configuration& = util::NoConfig());

    Trans(const Cache&, const FunctionSpace& gp, const FunctionSpace& sp,
          const eckit::Configuration& = util::NoConfig());
    Trans(const Cache&, const Grid&, int truncation, const eckit::Configuration& = util::NoConfig());
    Trans(const Cache&, const Grid&, const Domain&, int truncation, const eckit::Configuration& = util::NoConfig());

    void hash(eckit::Hash&) const;

    int truncation() const;
    size_t spectralCoefficients() const;
    const Grid& grid() const;
    const functionspace::Spectral& spectral() const;

    void dirtrans(const Field& gpfield, Field& spfield, const eckit::Configuration& = util::NoConfig()) const;

    void dirtrans(const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& = util::NoConfig()) const;

    void dirtrans_wind2vordiv(const Field& gpwind, Field& spvor, Field& spdiv,
                              const eckit::Configuration& = util::NoConfig()) const;

    void dirtrans_adj(const Field& spfield, Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;

    void dirtrans_adj(const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;

    void dirtrans_wind2vordiv_adj(const Field& spvor, const Field& spdiv, Field& gpwind,
                                  const eckit::Configuration& = util::NoConfig()) const;

    void invtrans(const Field& spfield, Field& gpfield, const eckit::Configuration& = util::NoConfig()) const;

    void invtrans(const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_grad(const Field& spfield, Field& gradfield, const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_grad(const FieldSet& spfields, FieldSet& gradfields,
                       const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_vordiv2wind(const Field& spvor, const Field& spdiv, Field& gpwind,
                              const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_adj(const Field& gpfield, Field& spfield, const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_adj(const FieldSet& gpfields, FieldSet& spfields,
                      const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_grad_adj(const Field& gradfield, Field& spfield,
                           const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_grad_adj(const FieldSet& gradfields, FieldSet& spfields,
                           const eckit::Configuration& = util::NoConfig()) const;

    void invtrans_vordiv2wind_adj(const Field& gpwind, Field& spvor, Field& spdiv,
                                  const eckit::Configuration& = util::NoConfig()) const;

    // -- IFS type fields --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
   * @brief invtrans
   * @param nb_scalar_fields
   * @param scalar_spectra
   * @param nb_vordiv_fields
   * @param vorticity_spectra
   * @param divergence_spectra
   * @param gp_fields
   */
    void invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                  const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                  const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief invtrans
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    void invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                  const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    void invtrans(const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                  double gp_fields[], const eckit::Configuration& = util::NoConfig()) const;


    /*!
   * @brief invtrans_adj
   * @param nb_scalar_fields
   * @param scalar_spectra
   * @param nb_vordiv_fields
   * @param vorticity_spectra
   * @param divergence_spectra
   * @param gp_fields
   */
    void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                      double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                      const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief invtrans_adj
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                      const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief Adjoint of Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    void invtrans_adj(const int nb_vordiv_fields, const double gp_fields[], double vorticity_spectra[],
                      double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief Direct transform of scalar fields
   */
    void dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                  const eckit::Configuration& = util::NoConfig()) const;

    /*!
   * @brief Direct transform of wind(U/V) to vorticity/divergence
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    void dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                  double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
