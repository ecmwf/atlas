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
#include <vector>


namespace atlas {
namespace grids { class Unstructured; }
class FieldGroup;
class Field;
namespace io {


class SimpleText {
public:

  /**
   * @brief Read SimpleText file into a Grid
   * @param file_path input file path
   * @param Grid data structure pointer to use
   * @return Grid data structure pointer
   */
  static grids::Unstructured* read(const std::string& file_path);

  /**
   * @brief Write Grid to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param grid Grid data structure
   */
  static void write(const std::string& file_path, const grids::Unstructured& grid);

  /**
   * @brief Write FieldGroup to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param fieldset FieldGroup data structure
   */
  static void write(const std::string& file_path, const FieldGroup& fieldset);

  /**
   * @brief Write Field to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param field Field data structure
   */
  static void write(const std::string& file_path, const Field& field);

  /**
   * @brief Write lan/lon and fields to SimpleText file (overwrites possibly existing file)
   * @param file_path output file path
   * @param nb_nodes number of points in unstructured grid
   * @param lon array (of length nb_nodes) pointer containing longitude (contiguous)
   * @param lat array (of length nb_nodes) pointer containing latitude (contiguous)
   * @param nb_fields number of fields (default none)
   * @param afvalues array (of length nb_fields) of arrays (of length nb_nodes), containing contiguous field values
   * @param afnames array (of length nb_fields) of field names
   */
  static void write(
      const std::string& file_path,
      const int& nb_nodes, const double*& lon, const double*& lat,
      const int& nb_fields=0, const double** afvalues=NULL , const char** afnames=NULL );

};


// C wrapper interfaces to C++ routines
extern "C"
{
  SimpleText*          atlas__simpletext__new ();
  void                 atlas__simpletext__delete        (SimpleText* This);
  grids::Unstructured* atlas__simpletext__read          (SimpleText* This, char* file_path);
  grids::Unstructured* atlas__read_simpletext           (char* file_path);
  void                 atlas__write_simpletext_fieldset (char* file_path, FieldGroup* fieldset);
  void                 atlas__write_simpletext_field    (char* file_path, Field* field);
}


} // namespace io
} // namespace atlas
#endif // SimpleText_h
