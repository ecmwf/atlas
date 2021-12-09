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
class Field;
}
}  // namespace atlas

namespace atlas {
namespace field {

// ------------------------------------------------------------------

/*!
 * \brief Field creator using array::ArrayShape parametrisation
 * \details
 * \code{.cpp}
 *    FieldImpl* field = Field::create(
 *         Config
 *           ("creator","ArraySpec")                // ArraySpec FieldCreator
 *           ("shape",array::make_shape(100,3))     // Rank 2 field with indexing [100][3]
 *           ("datatype",array::DataType::real64()) // Field internal data type
 *         );
 * \endcode
 */
class FieldCreatorArraySpec : public FieldCreator {
public:
    FieldCreatorArraySpec() {}
    FieldCreatorArraySpec(const eckit::Parametrisation&) {}
    virtual FieldImpl* createField(const eckit::Parametrisation&) const;
};

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
