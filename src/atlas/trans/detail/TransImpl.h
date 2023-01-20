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

#include "atlas/trans/Cache.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"

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

//-----------------------------------------------------------------------------

class TransImpl : public util::Object {
public:
    virtual std::string type() const { return "wrong value"; }

    virtual ~TransImpl() = 0;

    virtual int truncation() const = 0;

    virtual size_t nb_spectral_coefficients() const = 0;

    virtual size_t nb_spectral_coefficients_global() const = 0;

    virtual const Grid& grid() const = 0;

    virtual const functionspace::Spectral& spectral() const = 0;

    virtual void dirtrans(const Field& gpfield, Field& spfield,
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void dirtrans(const FieldSet& gpfields, FieldSet& spfields,
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void dirtrans_wind2vordiv(const Field& gpwind, Field& spvor, Field& spdiv,
                                      const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void dirtrans_adj(const Field& spfield, Field& gpfield,
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void dirtrans_adj(const FieldSet& spfields, FieldSet& gpfields,
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void dirtrans_wind2vordiv_adj(const Field& spvor, const Field& spdiv, Field& gpwind,
                                          const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans(const Field& spfield, Field& gpfield,
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans(const FieldSet& spfields, FieldSet& gpfields,
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_grad(const Field& spfield, Field& gradfield,
                               const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_grad(const FieldSet& spfields, FieldSet& gradfields,
                               const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_vordiv2wind(const Field& spvor, const Field& spdiv, Field& gpwind,
                                      const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_adj(const Field& gpfield, Field& spfield,
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_adj(const FieldSet& gpfields, FieldSet& spfields,
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_grad_adj(const Field& gradfield, Field& gpfield,
                                   const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_grad_adj(const FieldSet& gradfields, FieldSet& spfields,
                                   const eckit::Configuration& = util::NoConfig()) const = 0;

    virtual void invtrans_vordiv2wind_adj(const Field& gpwind, Field& spvor, Field& spdiv,
                                          const eckit::Configuration& = util::NoConfig()) const = 0;

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
    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                          const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief invtrans
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    virtual void invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    virtual void invtrans(const int nb_vordiv_fields, const double vorticity_spectra[],
                          const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief invtrans_adj
   * @param nb_scalar_fields
   * @param scalar_spectra
   * @param nb_vordiv_fields
   * @param vorticity_spectra
   * @param divergence_spectra
   * @param gp_fields
   */
    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                              double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief invtrans_adj
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    virtual void invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                              const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief Adjoint of Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    virtual void invtrans_adj(const int nb_vordiv_fields, const double wind_fields[], double vorticity_spectra[],
                              double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const = 0;


    /*!
   * @brief Direct transform of scalar fields
   */
    virtual void dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                          const eckit::Configuration& = util::NoConfig()) const = 0;

    /*!
   * @brief Direct transform of wind(U/V) to vorticity/divergence
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    virtual void dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                          double divergence_spectra[], const eckit::Configuration& = util::NoConfig()) const = 0;


    // deprecated
    virtual size_t spectralCoefficients() const { return nb_spectral_coefficients(); }

private:
    friend class TransInterface;
    // Only for TransIFS backend
    virtual int handle() const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
