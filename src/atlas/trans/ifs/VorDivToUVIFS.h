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

#include "atlas/trans/VorDivToUV.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class FunctionSpace;
}

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class VorDivToUVIFS : public trans::VorDivToUVImpl {
public:
    VorDivToUVIFS(const FunctionSpace&, const eckit::Configuration& = util::NoConfig());
    VorDivToUVIFS(int truncation, const eckit::Configuration& = util::NoConfig());

    virtual ~VorDivToUVIFS();

    virtual int truncation() const override { return truncation_; }

    // pure virtual interface

    // -- IFS style API --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
 * @brief Compute spectral wind (U/V) from spectral vorticity/divergence
 *
 * U = u*cos(lat)
 * V = v*cos(lat)
 *
 * @param nb_fields [in] Number of fields
 * @param vorticity [in] Spectral vorticity
 * @param divergence [in] Spectral divergence
 * @param U [out] Spectral wind U = u*cos(lat)
 * @param V [out] Spectral wind V = v*cos(lat)
 */
    virtual void execute(const int nb_coeff, const int nb_fields, const double vorticity[], const double divergence[],
                         double U[], double V[], const eckit::Configuration& = util::NoConfig()) const override;

private:
    int truncation_;
};

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
