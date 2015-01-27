/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef PointCloud_h
#define PointCloud_h


#include <string>
#include <vector>


// forward declarations
namespace eckit {
  class PathName;
}
namespace atlas {
  class FieldSet;
  class Field;
  class Grid;
  namespace grids { class Unstructured; }
}


namespace atlas {
namespace io {

/**
 * @brief PointCloud supports:
 * - reading atlas::grids::Unstructured
 * - writing of atlas::Grid
 * @warning only supports reading/writing doubles scalar fields
 */
 class PointCloud {
 public:

  /**
   * @brief Read PointCloud file into a Grid
   * @param path input file path
   * @param Grid data structure pointer to use
   * @return Grid data structure pointer
   */
  static grids::Unstructured* read(const eckit::PathName& path);

  /**
   * @brief Write Grid to PointCloud file (overwrites possibly existing file)
   * @param path output file path
   * @param grid Grid data structure
   */
  static void write(const eckit::PathName& path, const Grid& grid);

  /**
   * @brief Write FieldSet to PointCloud file (overwrites possibly existing file)
   * @param path output file path
   * @param fieldset FieldSet data structure
   */
  static void write(const eckit::PathName& path, const FieldSet& fieldset);

  /**
   * @brief Write lan/lon to PointCloud file (overwrites possibly existing file)
   * @note length of vectors lat and lon should be the same
   * @param path output file path
   * @param lon vector containing longitude values
   * @param lat vector containing latitude values
   */
  static void write(
      const eckit::PathName& path,
      const std::vector< Grid::Point >& pts );

  /**
   * @brief Write lan/lon and Field to PointCloud file (overwrites possibly existing file)
   * @param path output file path
   * @param field Field data structure
   */
  static void write(const eckit::PathName& path, const Field& field);

  /**
   * @brief Write lan/lon and fields to PointCloud file (overwrites possibly existing file)
   * @param path output file path
   * @param nb_pts number of points in unstructured grid
   * @param lon array (of length nb_pts) pointer containing longitude
   * @param lat array (of length nb_pts) pointer containing latitude
   * @param nb_fld number of fields (default none)
   * @param afvalues array (of length nb_fld) of arrays (of length nb_pts), containing field values
   * @param afnames array (of length nb_fld) of field names
   */
  static void write(
      const eckit::PathName& path,
      const int& nb_pts, const double* lon, const double* lat,
      const int& nb_fld=0, const double** afvalues=NULL , const char** afnames=NULL );

  /**
   * @brief Write lan/lon and fields to PointCloud file (overwrites possibly existing file)
   * @note length of vectors lat, lon and *vfvalues[] (each pointee) should be the same
   * @note length of vectors vfvalues and vfnames should be the same
   * @param path output file path
   * @param lon vector containing longitude values
   * @param lat vector containing latitude values
   * @param vfvalues vector of vector pointers, each pointing to field values vectors
   * @param vfnames vector of field names
   */
  static void write(
      const eckit::PathName& path,
      const std::vector< double >& lon, const std::vector< double >& lat,
      const std::vector< std::vector< double >* >& vfvalues=std::vector< std::vector< double >* >(),
      const std::vector< std::string >& vfnames=std::vector< std::string >() );

};


// C wrapper interfaces to C++ routines
extern "C"
{
  PointCloud*          atlas__pointcloud__new ();
  void                 atlas__pointcloud__delete        (PointCloud* This);
  grids::Unstructured* atlas__pointcloud__read          (PointCloud* This, char* file_path);
  grids::Unstructured* atlas__read_pointcloud           (char* file_path);
  void                 atlas__write_pointcloud_fieldset (char* file_path, FieldSet* fieldset);
  void                 atlas__write_pointcloud_field    (char* file_path, Field* field);
}


} // namespace io
} // namespace atlas
#endif // PointCloud_h
