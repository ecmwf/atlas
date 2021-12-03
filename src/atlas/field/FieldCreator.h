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

#ifndef atlas_field_FieldCreator_h
#define atlas_field_FieldCreator_h

#include <string>

#include "atlas/util/Object.h"

#include "atlas/field/Field.h"

namespace eckit {
class Parametrisation;
}

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace field {

//------------------------------------------------------------------------------------------------------

/*!
 * \brief Base class for creating new fields based on Parametrisation
 *
 * \details
 *  Example to create field[100][3] of default type double:
 * \code{.cpp}
 *    FieldImpl* field = Field::create(
 *         Config
 *           ("creator","ArraySpec")             // ArraySpec FieldCreator
 *           ("shape",array::make_shape(100,3))  // Rank 2 field with indexing [100][3]
 *         );
 * \endcode
 */
class FieldCreator : public util::Object {
public:
    FieldCreator();

    virtual ~FieldCreator();

    virtual FieldImpl* createField(const eckit::Parametrisation&) const = 0;
};

//------------------------------------------------------------------------------------------------------

class FieldCreatorFactory {
public:
    /*!
   * \brief build FieldCreator with factory key, and default options
   * \return FieldCreator
   */
    static FieldCreator* build(const std::string&);

    /*!
   * \brief build FieldCreator with options specified in parametrisation
   * \return mesh generator
   */
    static FieldCreator* build(const std::string&, const eckit::Parametrisation&);

    /*!
   * \brief list all registered field creators
   */
    static void list(std::ostream&);

private:
    std::string name_;
    virtual FieldCreator* make()                              = 0;
    virtual FieldCreator* make(const eckit::Parametrisation&) = 0;

protected:
    FieldCreatorFactory(const std::string&);
    virtual ~FieldCreatorFactory();
};

template <class T>
class FieldCreatorBuilder : public FieldCreatorFactory {
    virtual FieldCreator* make() { return new T(); }
    virtual FieldCreator* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    FieldCreatorBuilder(const std::string& name): FieldCreatorFactory(name) {}
};

//------------------------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas

#endif
