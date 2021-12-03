/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date June 2015

#pragma once

#include "atlas/field/FieldCreator.h"

namespace eckit {
class Parametrisation;
}
namespace atlas {
namespace field {
class FieldImpl;
}
}  // namespace atlas

namespace atlas {
namespace field {

// ------------------------------------------------------------------

/*!
 * \brief Field creator using IFS parametrisation
 * \details
 * Ideally this class should belong to IFS.
 * The only reference to IFS in Atlas should be here.
 * Example use:
 * \code{.cpp}
 *     FieldImpl* field = Field::create(
 *         Config
 *           ("creator","IFS")  // IFS FieldCreator
 *           ("ngptot",ngptot)  // Total number of grid points
 *           ("nproma",nproma)  // Grouping of grid points for vectorlength
 *           ("nlev",nlev)      // Number of levels
 *           ("nvar",nvar)      // Number of variables
 *           ("kind",8)         // Real kind in bytes
 *         );
 * \endcode
 */
class FieldCreatorIFS : public FieldCreator {
public:
    FieldCreatorIFS() {}
    FieldCreatorIFS(const eckit::Parametrisation&) {}
    virtual FieldImpl* createField(const eckit::Parametrisation&) const;
};

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
