/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Slavko Brdar
/// @date September 2024

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
 * \brief MultiField creator using datatype, shape, variable names as arguments
 * \details
 *     shape argument contains -1 at the position which gets filled with variable names
 * Example use:
 * \code{.cpp}
 *     MultiFieldImpl* multifield = MultiField::create(
 *         datatype,
 *         shape,
 *         var_names
 *         );
 * \endcode
 */
class MultiFieldCreatorArray : public MultiFieldCreator {
public:
    MultiFieldCreatorArray();
    MultiFieldCreatorArray(const eckit::Configuration& config);
    ~MultiFieldCreatorArray() override;
    MultiFieldImpl* create(const eckit::Configuration& config = util::Config()) const override;
    MultiFieldImpl* create(const array::DataType datatype, const array::ArrayShape& shape,
        const std::vector<std::string>& var_names) const override;
};

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
