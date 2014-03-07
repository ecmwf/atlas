// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#ifndef FieldSet_hpp
#define FieldSet_hpp

#include <vector>
#include <map>
#include <string>

namespace atlas {
class Field;

/// @brief Contains a list of field-pointers, no ownership
class FieldSet
{
public:
  FieldSet(const std::string& name="untitled");
  void add_field(Field& field);
  Field& field(const std::string& name);
  Field& field(size_t idx);
  std::vector<Field*>& fields() { return fields_; };
  size_t size() const { return fields_.size(); };
private:
  std::string name_;
  std::map< std::string, size_t > index_;
  std::vector< Field* > fields_;
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  FieldSet* atlas__FieldSet__new (char* name);
  void atlas__FieldSet__delete (FieldSet* This);
  void atlas__FieldSet__add_field (FieldSet* This, Field* field); 
  Field* atlas__FieldSet__field_by_name (FieldSet* This, char* name);
  int atlas__FieldSet__size (FieldSet* This);
  Field* atlas__FieldSet__field_by_idx (FieldSet* This, int idx);
  void atlas__FieldSet__fields (FieldSet* This, Field** &fields, int &nb_fields);
}
// ------------------------------------------------------------------

Field* atlas__FieldSet__field_by_name (FieldSet* This, char* name) {
  return &This->field( std::string(name) );
}

Field* atlas__FieldSet__field_by_idx (FieldSet* This, int idx) {
  return &This->field( idx );
}

} // namespace atlas

#endif // fieldset_hpp
