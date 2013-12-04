#ifndef fieldset_hpp
#define fieldset_hpp

#include <vector>
#include <map>
#include <string>

namespace ecmwf {
class Field;

/// @brief Contains a list of field-pointers, no ownership
class FieldSet
{
public:
  FieldSet(const std::string& name="untitled");
  void add_field(Field& field);
  Field& field(const std::string& name) { return *fields_[ index_.at(name) ]; }
  Field& field(size_t idx) { return *fields_[ idx ]; }
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
  FieldSet* ecmwf__FieldSet__new (char* name);
  void ecmwf__FieldSet__delete (FieldSet* This);
  void ecmwf__FieldSet__add_field (FieldSet* This, Field* field); 
  Field* ecmwf__FieldSet__field_by_name (FieldSet* This, char* name);
  int ecmwf__FieldSet__size (FieldSet* This);
  Field* ecmwf__FieldSet__field_by_idx (FieldSet* This, int idx);
  void ecmwf__FieldSet__fields (FieldSet* This, Field** &fields, int &nb_fields);
}
// ------------------------------------------------------------------

Field* ecmwf__FieldSet__field_by_name (FieldSet* This, char* name) {
  return &This->field( std::string(name) );
}

Field* ecmwf__FieldSet__field_by_idx (FieldSet* This, int idx) {
  return &This->field( idx );
}

} // namespace ecmwf

#endif // fieldset_hpp
