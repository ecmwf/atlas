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
/// @date Jan 2014

#ifndef atlas_FieldGroup_h
#define atlas_FieldGroup_h

#include <vector>
#include <map>
#include <string>

//------------------------------------------------------------------------------------------------------

namespace atlas {

class Field;

/// @brief Contains a list of field-pointers, no ownership

class FieldGroup {

public:

	FieldGroup(const std::string& name="untitled");

	const std::string& name() const { return name_; }

	void add_field(Field& field);
	bool has_field(const std::string& name) { return index_.count(name); }

	Field& field(const std::string& name);
	Field& field(size_t idx);

	std::vector<Field*>& fields() { return fields_; }
	size_t size() const { return fields_.size(); }

private:

	std::string name_;
	std::map< std::string, size_t > index_;
	std::vector< Field* > fields_;

};

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
	FieldGroup* atlas__FieldSet__new (char* name);
	void atlas__FieldSet__delete (FieldGroup* This);
	void atlas__FieldSet__add_field (FieldGroup* This, Field* field);
	int atlas__FieldSet__has_field (FieldGroup* This, char* name);
	Field* atlas__FieldSet__field_by_name (FieldGroup* This, char* name);
	int atlas__FieldSet__size (FieldGroup* This);
	Field* atlas__FieldSet__field_by_idx (FieldGroup* This, int idx);
	void atlas__FieldSet__fields (FieldGroup* This, Field** &fields, int &nb_fields);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_FieldGroup_h
