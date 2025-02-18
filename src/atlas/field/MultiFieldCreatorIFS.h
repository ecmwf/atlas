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

#include "atlas/field/MultiFieldCreator.h"

namespace eckit {
class Configuration;
}
namespace atlas {
namespace field {
class MultiFieldImpl;
}
}  // namespace atlas

namespace atlas {
namespace field {

// ------------------------------------------------------------------

/*!
 * \brief MultiField creator using IFS parametrisation
 * \details
 * Ideally this class should belong to IFS.
 * The only reference to IFS in Atlas::MultiField should be here.
 * Example use:
 * \code{.cpp}
 *     MultiFieldImpl* multifield = MultiField::create(
 *         Config
 *           ("creator","MultiFieldIFS")  // MultiFieldIFS FieldCreator
 *           ("ngptot",ngptot)  // Total number of grid points
 *           ("nproma",nproma)  // Grouping of grid points for vectorlength
 *           ("nlev",nlev)      // Number of levels
 *           ("nvar",nvar)      // Number of variables
 *           ("kind",8)         // Real kind in bytes
 *         );
 * \endcode
 */
class MultiFieldCreatorIFS : public MultiFieldCreator {
public:
    MultiFieldCreatorIFS();
    MultiFieldCreatorIFS(const eckit::Configuration& config);
    ~MultiFieldCreatorIFS() override;
    MultiFieldImpl* create(const eckit::Configuration& config = util::Config()) const override;
    MultiFieldImpl* create(const array::DataType datatype, const array::ArrayShape& shape,
        const std::vector<std::string>& var_names) const override;
};

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
