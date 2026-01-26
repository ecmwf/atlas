/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef atlas_field_MultiFieldCreator_h
#define atlas_field_MultiFieldCreator_h

#include <string>

#include "atlas/util/Object.h"

#include "atlas/field/MultiField.h"

namespace eckit {
class Configuration;
}

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace field {

//------------------------------------------------------------------------------------------------------

/*!
 * \brief Base class for creating new multifields based on Configuration
 *
 * \details
 *  Example to create field[100][3] of default type double:
 * \code{.cpp}
 *    FieldImpl* field = Field::create(
 *         Config
 *           ("creator","ArraySpec")             // ArraySpec MultiFieldCreator
 *           ("shape",array::make_shape(100,3))  // Rank 2 field with indexing [100][3]
 *         );
 * \endcode
 */
class MultiFieldCreator : public util::Object {
public:
    MultiFieldCreator();
    MultiFieldCreator(const eckit::Configuration& config);

    virtual ~MultiFieldCreator();

    virtual MultiFieldImpl* create(const eckit::Configuration& config = util::Config()) const = 0;
    virtual MultiFieldImpl* create(const array::DataType datatype, const array::ArrayShape& shape,
        const std::vector<std::string>& var_names) const = 0;
};

//------------------------------------------------------------------------------------------------------

class MultiFieldCreatorFactory : public util::Factory<MultiFieldCreatorFactory> {
public:
    static std::string className() { return "MultiFieldCreatorFactory"; }

    /*!
   * \brief build MultiFieldCreator with options specified in parametrisation
   * \return mesh generator
   */
    static MultiFieldCreator* build(const std::string&, const eckit::Configuration& = util::Config());

    using Factory::Factory;

private:
    virtual MultiFieldCreator* make()                            = 0;
    virtual MultiFieldCreator* make(const eckit::Configuration&) = 0;
};

template <class T>
class MultiFieldCreatorBuilder : public MultiFieldCreatorFactory {
    virtual MultiFieldCreator* make() { return new T(); }
    virtual MultiFieldCreator* make(const eckit::Configuration& config) { return new T(config); }

public:
    using MultiFieldCreatorFactory::MultiFieldCreatorFactory;
};

//------------------------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas

#endif
