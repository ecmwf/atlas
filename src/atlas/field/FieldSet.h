/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Pedro Maciel
/// @date Jan 2015

#pragma once

#include <vector>
#include <iterator>


#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/field/Field.h"

namespace atlas {
namespace field {

class FieldSet;

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSetImpl : public eckit::Owned {

public: // types

  typedef std::vector<Field>::iterator       iterator;
  typedef std::vector<Field>::const_iterator const_iterator;

public: // methods

  /// Constructs an empty FieldSet
  FieldSetImpl(const std::string& name = "untitled");

  size_t size() const { return  fields_.size(); }
  bool empty()  const { return !fields_.size(); }

  void clear();

  const std::string& name() const { return name_; }
        std::string& name()       { return name_; }

  const Field& operator[](const size_t &i) const { return field(i); }
        Field& operator[](const size_t &i)       { return field(i); }

  const Field& operator[](const std::string &name) const { return field(name); }
        Field& operator[](const std::string &name)       { return field(name); }

  const Field& field(const size_t& i) const { ASSERT(i<size()); return fields_[i]; }
        Field& field(const size_t& i)       { ASSERT(i<size()); return fields_[i]; }

  std::vector< std::string > field_names() const;

  Field add(const Field&);

  bool has_field(const std::string& name) const;

  Field& field(const std::string& name) const;

  iterator begin() { return fields_.begin(); }
  iterator end()   { return fields_.end(); }
  const_iterator begin()  const { return fields_.begin(); }
  const_iterator end()    const { return fields_.end(); }
  const_iterator cbegin() const { return fields_.begin(); }
  const_iterator cend()   const { return fields_.end(); }

protected: // data

  std::vector< Field >             fields_;  ///< field storage
  std::string                      name_;    ///< internal name
  std::map< std::string, size_t >  index_;   ///< name-to-index map, to refer fields by name
};

class FieldSet;

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSet : public eckit::Owned {

public: // types

  typedef std::vector<Field>::iterator       iterator;
  typedef std::vector<Field>::const_iterator const_iterator;

public: // methods

  FieldSet( const std::string& name = "untitled" );
  FieldSet( const FieldSetImpl* );
  FieldSet( const FieldSet& );

  size_t size() const { return  fieldset_->size(); }
  bool empty()  const { return  fieldset_->empty(); }

  void clear() { fieldset_->clear(); }

  const std::string& name() const { return fieldset_->name(); }
        std::string& name()       { return fieldset_->name(); }

  const Field& operator[](const size_t &i) const { return fieldset_->operator[](i); }
        Field& operator[](const size_t &i)       { return fieldset_->operator[](i); }

  const Field& operator[](const std::string &name) const { return fieldset_->operator[](name); }
        Field& operator[](const std::string &name)       { return fieldset_->operator[](name); }

  const Field& field(const size_t& i) const { return fieldset_->field(i); }
        Field& field(const size_t& i)       { return fieldset_->field(i); }

  std::vector< std::string > field_names() const { return fieldset_->field_names(); }

  Field add( const Field& field ) { return fieldset_->add(field); }

  bool has_field(const std::string& name) const { return fieldset_->has_field(name); }

  Field& field(const std::string& name) const { return fieldset_->field(name); }

  iterator begin() { return fieldset_->begin(); }
  iterator end()   { return fieldset_->end(); }
  const_iterator begin()  const { return fieldset_->begin(); }
  const_iterator end()    const { return fieldset_->end(); }
  const_iterator cbegin() const { return fieldset_->begin(); }
  const_iterator cend()   const { return fieldset_->end(); }

private: // data

  eckit::SharedPtr<FieldSetImpl> fieldset_;
};


// C wrapper interfaces to C++ routines
extern "C"
{
  FieldSetImpl* atlas__FieldSet__new       (char* name);
  void      atlas__FieldSet__delete        (FieldSetImpl* This);
  void      atlas__FieldSet__add_field     (FieldSetImpl* This, FieldImpl* field);
  int       atlas__FieldSet__has_field     (const FieldSetImpl* This, char* name);
  size_t    atlas__FieldSet__size          (const FieldSetImpl* This);
  FieldImpl*    atlas__FieldSet__field_by_name (FieldSetImpl* This, char* name);
  FieldImpl*    atlas__FieldSet__field_by_idx  (FieldSetImpl* This, size_t idx);
}

} // namespace field
} // namespace atlas
