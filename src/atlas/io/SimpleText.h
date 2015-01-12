/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SimpleText_h
#define SimpleText_h

#include <string>


namespace atlas {
class Mesh;
class FieldGroup;
class Field;
namespace io {


class SimpleText {
public:

  /**
   * @brief Read SimpleText file
   * @param file_path input file path
   * @param mesh data structure pointer to reuse, allocates new data structure if NULL (default)
   * @return mesh data structure pointer
   */
  static Mesh* read(const std::string& file_path, Mesh* mesh=NULL);

  /**
   * @brief Write SimpleText file
   * @param file_path output file path
   * @param mesh mesh data structure
   */
  void write(const std::string& file_path, const Mesh& mesh) const;

  /**
   * @brief Write FieldGroup to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param fieldset FieldGroup data structure
   */
  void write(const std::string& file_path, const FieldGroup& fieldset) const;

  /**
   * @brief Write Field to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param field Field data structure
   */
  void write(const std::string& file_path, const Field& field) const;

};


// C wrapper interfaces to C++ routines
extern "C"
{
  SimpleText* atlas__simpletext__new ();
  void        atlas__simpletext__delete        (SimpleText* This);
  Mesh*       atlas__simpletext__read          (SimpleText* This, char* file_path);
  void        atlas__simpletext__write         (SimpleText* This, char* file_path, Mesh* mesh);
  Mesh*       atlas__read_simpletext           (char* file_path);
  void        atlas__write_simpletext_mesh     (char* file_path, Mesh* mesh);
  void        atlas__write_simpletext_fieldset (char* file_path, FieldGroup* fieldset);
  void        atlas__write_simpletext_field    (char* file_path, Field* field);
}


} // namespace io
} // namespace atlas

#endif // SimpleText_h
