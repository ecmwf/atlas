/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/util/Point.h"

// forward declarations
namespace eckit {
class PathName;
}

namespace atlas {
class Field;
class FieldSet;
}  // namespace atlas

namespace atlas {
namespace grid {
class Grid;
class Unstructured;
}  // namespace grid
}  // namespace atlas

namespace atlas {
class Mesh;
}
namespace atlas {
namespace functionspace {
class NodeColumns;
}
}  // namespace atlas

namespace atlas {
namespace output {
namespace detail {

/**
 * @brief PointCloudIO supports:
 * - reading Mesh
 * - writing of Mesh
 * @warning only supports reading/writing doubles scalar fields
 */
class PointCloudIO {
public:
    /**
  * @brief Read PointCloudIO file into a Mesh
  * @param path input file path
  * @param Grid data structure pointer to use
  * @return Grid data structure pointer
  */
    static Mesh read(const eckit::PathName& path);

    /**
  * @brief Read PointCloudIO file into a Mesh
  * @param path input file path
  * @param vfnames names of fields to read
  * @return Mesh pointer
  */
    static Mesh read(const eckit::PathName& path, std::vector<std::string>& vfnames);

    /**
 * @brief Write Grid to PointCloudIO file (overwrites possibly existing file)
 * @param path output file path
 * @param grid Grid data structure
 */
    static void write(const eckit::PathName& path, const Mesh& mesh);

    /**
 * @brief Write FieldSet to PointCloudIO file (overwrites possibly existing
 * file)
 * @param path output file path
 * @param fieldset FieldSet data structure
 */
    static void write(const eckit::PathName& path, const FieldSet& fieldset,
                      const functionspace::NodeColumns& function_space);

    /**
 * @brief Write lan/lon to PointCloudIO file (overwrites possibly existing file)
 * @note length of vectors lat and lon should be the same
 * @param path output file path
 * @param lon vector containing longitude values
 * @param lat vector containing latitude values
 */
    static void write(const eckit::PathName& path, const std::vector<PointLonLat>& pts);

    /**
 * @brief Write lan/lon and fields to PointCloudIO file (overwrites possibly
 * existing file)
 * @param path output file path
 * @param nb_pts number of points in unstructured grid
 * @param lon array (of length nb_pts) pointer containing longitude
 * @param lat array (of length nb_pts) pointer containing latitude
 * @param nb_fld number of fields (default none)
 * @param afvalues array (of length nb_fld) of arrays (of length nb_pts),
 * containing field values
 * @param afnames array (of length nb_fld) of field names
 */
    static void write(const eckit::PathName& path, const int& nb_pts, const double* lon, const double* lat,
                      const int& nb_fld = 0, const double** afvalues = NULL, const char** afnames = NULL);

    /**
 * @brief Write lan/lon and fields to PointCloudIO file (overwrites possibly
 * existing file)
 * @note length of vectors lat, lon and *vfvalues[] (each pointee) should be the
 * same
 * @note length of vectors vfvalues and vfnames should be the same
 * @param path output file path
 * @param lon vector containing longitude values
 * @param lat vector containing latitude values
 * @param vfvalues vector of vector pointers, each pointing to field values
 * vectors
 * @param vfnames vector of field names
 */
    static void write(const eckit::PathName& path, const std::vector<double>& lon, const std::vector<double>& lat,
                      const std::vector<std::vector<double>*>& vfvalues = std::vector<std::vector<double>*>(),
                      const std::vector<std::string>& vfnames           = std::vector<std::string>());
};

// C wrapper interfaces to C++ routines
extern "C" {
// PointCloudIO* atlas__pointcloud__new ();
// void atlas__pointcloud__delete(PointCloudIO* This);
// mesh::Mesh* atlas__pointcloud__read(PointCloudIO* This, char* file_path);
// mesh::Mesh* atlas__read_pointcloud(char* file_path);
// void atlas__write_pointcloud_fieldset(char* file_path, const
// field::FieldImpl* fieldset, const functionspace::FunctionSpaceImpl*
// functionspace);
// void atlas__write_pointcloud_field(char* file_path, const field::FieldImpl*
// field, const functionspace::FunctionSpaceImpl* functionspace);
}

}  // namespace detail
}  // namespace output
}  // namespace atlas
