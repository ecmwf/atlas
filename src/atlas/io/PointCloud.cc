/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Log.h"

#include "atlas/Field.h"
#include "atlas/FieldSet.h"
#include "atlas/FunctionSpace.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/util/ArrayView.h"

#include "atlas/io/PointCloud.h"


namespace atlas {
namespace io {


// ------------------------------------------------------------------


namespace {


std::string sanitize_field_name(const std::string& s)
{
  // replace non-printable characters, then trim right & left
  std::string r(s);
  std::replace_if(r.begin(), r.end(), ::isspace, '_');
  r.erase(r.find_last_not_of('_')+1);
  r.erase(0, r.find_first_not_of('_'));
  if (!r.length())
    r = "_";
  return r;
}


} // end anonymous namespace


// ------------------------------------------------------------------


grids::Unstructured* PointCloud::read(const eckit::PathName& path, std::vector<std::string>& vfnames )
{
  const std::string msg("PointCloud::read: ");

  // read file into a vector of points and a field matrix
  // @warning: several copy operations here

  vfnames.clear();

  std::vector< std::vector< double > > vfvalues;  // field values
  std::vector< Grid::Point >* pts = 0;            // (lon,lat) points, to build unstructured grid

  {
    std::string line;
    std::istringstream iss;
    std::ostringstream oss;
    size_t nb_pts;      // # points/nodes
    size_t nb_columns;  // # data columns (lon,lat,field1,field2,...)
    size_t nb_fld;      // # fields (nb_fld := nb_columns-2)

    // open file and read all of header & data
    std::ifstream f(path.asString().c_str());
    if (!f.is_open())
      throw eckit::CantOpenFile(path.asString());

    // header, part 1:
    // determine number of rows/columns
    // (read all of line, look for "PointCloud" signature, nb_pts, nb_columns, ...)
    std::getline(f,line);
    iss.str(line);
    iss >> line >> nb_pts >> nb_columns;
    if (line!="PointCloud")
      throw eckit::BadParameter(msg+"beggining of file not found (expected: PointCloud)");
    if (nb_pts==0)
      throw eckit::BadValue(msg+"invalid number of points (failed: nb_pts>0)");
    if (nb_columns<2)
      throw eckit::BadValue(msg+"invalid number of columns (failed: nb_columns>=2)");


    // header, part 2:
    // determine columns' labels
    // (check end of first line for possible column labels, starting from defaults)

    vfnames.resize(nb_columns);
    for (size_t j=0; j<nb_columns; ++j) {
      oss.str("column_");
      oss << (j+1);
      vfnames[j] = (iss && iss>>line)? sanitize_field_name(line) : oss.str();
    }


    // (preallocate data, and define fields without the first two columns (lon,lat) because they are treated differently)
    vfnames.erase(vfnames.begin(),vfnames.begin()+2);
    nb_fld = nb_columns-2;  // always >= 0, considering previous check
    if (nb_fld)
      vfvalues.assign(nb_fld,std::vector< double >(nb_pts,0.));


    // data
    // (grid is created after reading coordinates because the internal bounding box needs them)

    pts = new std::vector< Grid::Point >(nb_pts);

    size_t i,j;  // (index for node/row and field/column, out of scope to check at end of loops)
    for (i=0; f && i<nb_pts; ++i)
    {
      std::getline(f,line);
      iss.clear();
      iss.str(line);
      double lon, lat;

      //NOTE always expects (lon,lat) order, maybe make it configurable?
      iss >> lon >> lat;
      for (j=0; iss && j<nb_fld; ++j)
        iss >> vfvalues[j][i];
      if (j<nb_fld) {
        oss << "invalid number of fields in data section, on line " << (i+1) << ", read " << j << " fields, expected " << nb_fld << ".";
        throw eckit::BadValue(msg+oss.str());
      }

      pts->operator[](i).assign(lon,lat);
    }
    if (i<nb_pts) {
      oss << "invalid number of lines in data section, read " << (i) << " lines, expected " << nb_pts << ".";
      throw eckit::BadValue(msg+oss.str());
    }

    f.close();
  }


  // build unstructured grid, to return
  // (copies read-in data section into scalar fields)

  grids::Unstructured* grid = new grids::Unstructured(pts);
  ASSERT(grid);

  Mesh& m = grid->mesh();
  for (size_t j=0; j<vfvalues.size(); ++j) {
    Field& field = m.function_space("nodes").create_field< double >(vfnames[j],1);
    ArrayView< double, 1 > fdata(field);
    for (size_t i=0; i<grid->npts(); ++i)
      fdata(i) = vfvalues[j][i];
  }

  return grid;
}


grids::Unstructured*PointCloud::read(const eckit::PathName& path)
{
  std::vector<std::string> vfnames;
  return read(path,vfnames);
}

void PointCloud::write(const eckit::PathName& path, const Grid& grid)
{
  const std::string msg("PointCloud::write: ");

  // operate in grid's mesh function space, creating transversing data structures
  // @warning: several copy operations here

  FunctionSpace& nodes = grid.mesh().function_space("nodes");

  ArrayView< double, 2 > lonlat(nodes.field("lonlat"));
  if (!lonlat.size())
    throw eckit::BadParameter(msg+"invalid number of points (failed: nb_pts>0)");

  // get the fields (sanitized) names and values
  // (bypasses fields ("lonlat"|"coordinates") as shape(1)!=1)
  std::vector< std::string > vfnames;
  std::vector< ArrayView< double, 1 > > vfvalues;
  for (int i=0; i<nodes.nb_fields(); ++i)
  {
    const Field& field = nodes.field(i);
    if ( field.shape(0)==lonlat.size() &&
         field.shape(1)==1 &&
         field.name()!="glb_idx" )  // FIXME: no support for non-int types!
    {
      vfnames.push_back(sanitize_field_name(field.name()));
      vfvalues.push_back(ArrayView< double, 1 >(field));
    }
  }

  std::ofstream f( path.asString().c_str() );
  if (!f.is_open())
    throw eckit::CantOpenFile(path.asString());

  const size_t Npts = lonlat.size();
  const size_t Nfld = vfvalues.size();

  // header
  f << "PointCloud\t" << Npts << '\t' << (2+Nfld) << "\tlon\tlat";
  for (size_t j=0; j<Nfld; ++j)
    f << '\t' << vfnames[j];
  f << '\n';

  // data
  for (size_t i=0; i<Npts; ++i) {
    f << lonlat(i,0) << '\t' << lonlat(i,1);
    for (size_t j=0; j<Nfld; ++j)
      f << '\t' << vfvalues[j](i);
    f << '\n';
  }

  f.close();
}


void PointCloud::write(const eckit::PathName& path, const FieldSet& fieldset)
{
  const std::string msg("PointCloud::write: ");

  // operate in field sets with same grid and consistent size(s), creating transversing data structures
  // @warning: several copy operations here

  for (size_t i=1; i<fieldset.size(); ++i)
    if (!fieldset[0].grid().same( fieldset[i].grid() ))
      throw eckit::BadParameter(msg+"fields must be described in the same grid (fieldset.field(0).grid() == fieldset.field(*).grid())");

  const Mesh& m(fieldset[0].grid().mesh());
  ArrayView< double, 2 > lonlat(m.function_space("nodes").field("lonlat"));
  if (!lonlat.size())
    throw eckit::BadParameter(msg+"invalid number of points (failed: nb_pts>0)");

  // get the fields (sanitized) names and values
  // (bypasses fields ("lonlat"|"coordinates") as shape(1)!=1)
  std::vector< std::string > vfnames;
  std::vector< ArrayView< double, 1 > > vfvalues;
  for (size_t i=0; i<fieldset.size(); ++i)
  {
    const Field& field = fieldset[i];
    if ( field.shape(0)==lonlat.size() &&
         field.shape(1)==1 &&
         field.name()!="glb_idx" )  // FIXME: no support for non-int types!
    {
      vfnames.push_back(sanitize_field_name(field.name()));
      vfvalues.push_back(ArrayView< double, 1 >(field));
    }
  }

  std::ofstream f(path.asString().c_str());
  if (!f.is_open())
    throw eckit::CantOpenFile(path.asString());
  const size_t
      Npts = lonlat.size(),
      Nfld = vfvalues.size();

  // header
  f << "PointCloud\t" << Npts << '\t' << (2+Nfld) << "\tlon\tlat";
  for (size_t j=0; j<Nfld; ++j)
    f << '\t' << vfnames[j];
  f << '\n';

  f.precision( std::numeric_limits< double >::digits10 );

  // data
  for (size_t i=0; i<Npts; ++i) {
    f << lonlat(i,0) << '\t' << lonlat(i,1);
    for (size_t j=0; j<Nfld; ++j)
      f << '\t' << vfvalues[j](i);
    f << '\n';
  }

  f.close();
}


void PointCloud::write(
    const eckit::PathName& path,
    const std::vector< Grid::Point >& pts )
{
  std::ofstream f(path.asString().c_str());
  if (!f.is_open())
    throw eckit::CantOpenFile(path.asString());

  // header
  f << "PointCloud\t" << pts.size() << '\t' << 2 << "\tlon\tlat\n";

  // data
  for (size_t i=0; i<pts.size(); ++i)
    f << pts[i].lon() << '\t' << pts[i].lat() << '\n';

  f.close();
}


void PointCloud::write(const eckit::PathName& path, const Field& field)
{
  FieldSet fieldset;
  fieldset.add_field(const_cast< Field& >(field));
  write(path,fieldset);
}


void PointCloud::write(
    const eckit::PathName& path,
    const std::vector< double >& lon, const std::vector< double >& lat,
    const std::vector< std::vector< double >* >& vfvalues,
    const std::vector< std::string >& vfnames )
{
  const std::string msg("PointCloud::write: ");
  const size_t
      Npts (lon.size()),
      Nfld (vfvalues.size());
  if (Npts!=lat.size())
    throw eckit::BadParameter(msg+"number of points inconsistent (failed: #lon == #lat)");
  if (Nfld!=vfnames.size())
    throw eckit::BadParameter(msg+"number of fields inconsistent (failed: #vfvalues == #vfnames)");
  for (size_t j=0; j<Nfld; ++j)
    if (Npts!=vfvalues[j]->size())
      throw eckit::BadParameter(msg+"number of points inconsistent (failed: #lon == #lat == #*vfvalues[])");

  std::ofstream f(path.asString().c_str());
  if (!f.is_open())
    throw eckit::CantOpenFile(path.asString());

  // header
  f << "PointCloud\t" << Npts << '\t' << (2+Nfld) << "\tlon\tlat";
  for (size_t j=0; j<Nfld; ++j)
    f << '\t' << sanitize_field_name(vfnames[j]);
  f << '\n';

  // data
  for (size_t i=0; i<Npts; ++i) {
    f << lon[i] << '\t' << lat[i];
    for (size_t j=0; j<Nfld; ++j)
      f << '\t' << vfvalues[j]->operator[](i);
    f << '\n';
  }

  f.close();
}


void PointCloud::write(
    const eckit::PathName& path,
    const int& nb_pts, const double* lon, const double* lat,
    const int& nb_fld, const double** afvalues, const char** afnames )
{
  const std::string msg("PointCloud::write: ");

  const size_t
      Npts (nb_pts>0? nb_pts : 0),
      Nfld (nb_fld>0 && afvalues && afnames? nb_fld : 0);
  if (!Npts)
    throw eckit::BadParameter(msg+"invalid number of points (nb_nodes)");
  if (!lon)
    throw eckit::BadParameter(msg+"invalid array describing longitude (lon)");
  if (!lat)
    throw eckit::BadParameter(msg+"invalid array describing latitude (lat)");

  std::ofstream f(path.asString().c_str());
  if (!f.is_open())
    throw eckit::CantOpenFile(path.asString());

  // header
  f << "PointCloud\t" << Npts << '\t' << (2+Nfld) << "\tlon\tlat";
  for (size_t j=0; j<Nfld; ++j)
    f << '\t' << sanitize_field_name(afnames[j]);
  f << '\n';

  // data
  for (size_t i=0; i<Npts; ++i) {
    f << lon[i] << '\t' << lat[i];
    for (size_t j=0; j<Nfld; ++j)
      f << '\t' << afvalues[j][i];
    f << '\n';
  }

  f.close();
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines


PointCloud* atlas__pointcloud__new()
{ return new PointCloud(); }


void atlas__pointcloud__delete (PointCloud* This)
{ delete This; }


grids::Unstructured* atlas__pointcloud__read (PointCloud* This, char* file_path)
{ return This->read(file_path); }


grids::Unstructured* atlas__read_pointcloud (char* file_path)
{ return PointCloud::read(file_path); }


void atlas__write_pointcloud_fieldset (char* file_path, FieldSet* fieldset)
{ PointCloud::write(file_path, *fieldset); }


void atlas__write_pointcloud_field (char* file_path, Field* field)
{ PointCloud::write(file_path, *field); }


// ------------------------------------------------------------------


} // namespace io
} // namespace atlas
