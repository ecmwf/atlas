/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date June 2015

#ifndef atlas_field_FieldTCreator_h
#define atlas_field_FieldTCreator_h

#include <iosfwd>
#include <string>
#include "eckit/memory/Owned.h"
#include "atlas/field/FieldTCreator.h"
#include "atlas/util/ArrayUtil.h"

namespace eckit { class Parametrisation; }
namespace atlas { class Field; }

namespace atlas {
namespace field {

// ------------------------------------------------------------------

/*!
 * \brief FieldT<DATA_TYPE> creator using ArrayShape and parametrisation
 * \code{.cpp}
 *    Field* field = Field::create(
 *         Field::Parameters
 *           ("creator","ArraySpec")     // ArraySpec FieldCreator
 *           ("shape",make_shape(100,3))  // Rank 2 field with indexing [100][3]
 *           ("datatype",DataType::real64()) // Field internal data type
 *         );
 * \endcode
 */

class FieldTCreator: public eckit::Owned
{
public:
  FieldTCreator();
  virtual ~FieldTCreator();
  virtual Field* create_field( const ArrayShape&, const eckit::Parametrisation& ) const = 0;
};

// ------------------------------------------------------------------

class FieldTCreatorFactory {
  public:
    /*!
     * \brief build FieldCreator with factory key, and default options
     * \return FieldCreator
     */
    static FieldTCreator* build(const std::string&);

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);

  private:
    std::string name_;
    virtual FieldTCreator* make() = 0 ;

  protected:

    FieldTCreatorFactory(const std::string&);
    virtual ~FieldTCreatorFactory();
};

// ------------------------------------------------------------------

template<class T>
class FieldTCreatorBuilder : public FieldTCreatorFactory {
  virtual FieldTCreator* make() {
      return new T();
  }
  public:
    FieldTCreatorBuilder(const std::string& name) : FieldTCreatorFactory(name) {}
};

// ------------------------------------------------------------------

template <typename DATA_TYPE>
class FieldTCreatorT: public FieldTCreator
{
public:
  FieldTCreatorT() {}
  virtual ~FieldTCreatorT() {}
  virtual Field* create_field( const ArrayShape&, const eckit::Parametrisation& ) const;
};

// ------------------------------------------------------------------

} // namespace field
} // namespace atlas

#endif

