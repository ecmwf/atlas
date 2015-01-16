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
#include <algorithm>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/log/Log.h"
#include "atlas/Field.h"
#include "atlas/FieldGroup.h"
#include "atlas/FunctionSpace.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/util/ArrayView.h"
#include "atlas/io/SimpleText.h"


//using namespace eckit;
using eckit::LocalPathName;
using eckit::Log;


//NOTE: this class only supports reading/writing doubles


namespace atlas {
namespace io {


// ------------------------------------------------------------------


namespace {


template< typename DATA_TYPE >
void write_field_nodes(const Field& field, std::ostream& out)
{
#if 0
  ArrayView< DATA_TYPE > data(field);
#endif
}


void sanitize_field_name(std::string& s)
{
  // replace non-printable characters, then trim right & left
  std::replace_if(s.begin(), s.end(), ::isspace, '_');
  s.erase(s.find_last_not_of('_')+1);
  s.erase(0, s.find_first_not_of('_'));
  if (!s.length())
    s = "_";
}


} // end anonymous namespace


// ------------------------------------------------------------------


grids::Unstructured* SimpleText::read(const std::string& file_path)
{
  std::vector< Grid::Point >* pts = 0;
  std::vector< std::string > fnames;
  std::vector< std::vector< double > > fvalues;

  std::ostringstream oss;
  oss << "SimpleText::read: ";

  {
    LocalPathName path(file_path);
    std::ifstream file(path.c_str());
    if (!file.is_open())
      throw eckit::CantOpenFile(path);


    // read header section
    // (validates nb_pts, nb_columns, nb_fields, and fnames)

    std::string discard;
    size_t
        nb_pts,
        nb_columns,
        i,
        j;

    file >> nb_pts >> nb_columns;
    std::getline(file,discard);
    const size_t nb_fld = (nb_columns>2? nb_columns-2 : 0);

    Log::info() << "SimpleText::read: nb_pts/nb_columns: " << nb_pts << " / " << nb_columns << std::endl;
    if (nb_pts==0)
      throw eckit::BadValue("SimpleText::read: invalid number of points (nb_pts==0)");
    if (nb_columns<2)
      throw eckit::BadValue("SimpleText::read: invalid number of columns (nb_columns<2)");

    file >> discard >> discard;  // (skip "lat lon")
    if (nb_fld)
      fnames.resize(nb_fld);
    for (j=0; file && j<nb_fld; ++j)
      file >> fnames[j];
    std::getline(file,discard);// (discard until end of line)

    if (j!=nb_fld) {
      oss << "invalid number of fields in header section, "
             "read " << j << " fields, expected " << nb_fld << ".";
      throw eckit::BadValue(oss.str());
    }


    // read data section
    // (grid is created after reading coordinates because the internal bounding box needs them)

    pts = new std::vector< Grid::Point >(nb_pts);
    if (nb_fld)
      fvalues.assign(nb_fld,std::vector< double >(nb_pts,0.));

    for (i=0; file && i<nb_pts; ++i) {
      double lon, lat;

      file >> lon >> lat;
      pts->operator[](i).assign(lon,lat);

      for (j=0; file && j<nb_fld; ++j)
        file >> fvalues[j][i];
      if (j<nb_fld) {
        oss << "invalid number of fields in data section, "
               "on line " << i << ", read " << j << " fields, expected " << nb_fld << ".";
        throw eckit::BadValue(oss.str());
      }

      std::getline(file,discard);  // (discard until end of line)
    }
    if (i<nb_pts) {
      oss << "invalid number of lines in data section, "
             "read " << i << " lines, expected " << nb_pts << ".";
      throw eckit::BadValue(oss.str());
    }


    file.close();
  }


  // build unstructured grid, to return
  // (copies read-in data section into scalar fields)

  grids::Unstructured* grid = new grids::Unstructured(pts);
  ASSERT(grid);
  Mesh& m = grid->mesh();  
  for (size_t j=0; j<fvalues.size(); ++j) {
    Field& f = m.function_space("nodes").create_field< double >(fnames[j],1);
    ArrayView< double, 1 > field(f);
    for (size_t i=0; i<grid->npts(); ++i)
      field(i) = fvalues[j][i];
  }

  return grid;
}


void SimpleText::write(const std::string &file_path, const grids::Unstructured &grid)
{
  NOTIMP;
}


void SimpleText::write(const std::string& file_path, const FieldGroup &fieldset)
{
  LocalPathName path(file_path);
  Log::info() << "SimpleText::write: FieldGroup " << fieldset.name() << " to file " << path << std::endl;
  std::ofstream file(path.c_str());

  for (int i=0; i<fieldset.size(); ++i) {
    const Field& field = fieldset.field(i);
    if (field.function_space().metadata().has< int >("type") &&
        field.function_space().metadata().get< int >("type") == Entity::NODES) {
      field.data_type()=="int32"?  write_field_nodes< int    >(field,file) :
      field.data_type()=="int64"?  write_field_nodes< long   >(field,file) :
      field.data_type()=="real32"? write_field_nodes< float  >(field,file) :
      field.data_type()=="real64"? write_field_nodes< double >(field,file) :
                                   void();
    }
    else {
      throw eckit::Exception("SimpleText::write: function_space "+field.function_space().name()+" has no type.. ?");
    }
    file << std::flush;
  }

  file.close();
}


void SimpleText::write(const std::string& file_path, const Field &field)
{
  LocalPathName path(file_path);
  Log::info() << "SimpleText::write: Field " << field.name() << " to file " << path << std::endl;
  std::ofstream file(path.c_str());

  if (field.function_space().metadata().has< int >("type") &&
      field.function_space().metadata().get< int >("type") == Entity::NODES) {
    field.data_type()=="int32"?  write_field_nodes< int    >(field,file) :
    field.data_type()=="int64"?  write_field_nodes< long   >(field,file) :
    field.data_type()=="real32"? write_field_nodes< float  >(field,file) :
    field.data_type()=="real64"? write_field_nodes< double >(field,file) :
                                 void();
  }
  else {
    throw eckit::Exception("SimpleText::write: function_space "+field.function_space().name()+" has no type.?");
  }

  file << std::flush;
  file.close();
}


void SimpleText::write(
    const std::string& file_path,
    const int& nb_pts, const double* lon, const double* lat,
    const int& nb_fld, const double** afvalues, const char** afnames )
{
  LocalPathName path(file_path);

  const size_t
      Npts (nb_pts>0? nb_pts : 0),
      Nfld (nb_fld>0? nb_fld : 0);
  if (!Npts)
    throw eckit::BadParameter("SimpleText::write: invalid number of points (nb_nodes)");
  if (!lon)
    throw eckit::BadParameter("SimpleText::write: invalid array describing longitude (lon)");
  if (!lat)
    throw eckit::BadParameter("SimpleText::write: invalid array describing latitude (lat)");

  std::ofstream file(path.c_str());
  if (!file.is_open())
    throw eckit::CantOpenFile(path);

  // header
  file << Npts << '\t' << (2+Nfld) << "\n"
          "lon\tlat";
  for (size_t j=0; j<Nfld; ++j) {
    std::string fname(afnames[j]);
    sanitize_field_name(fname);
    file << '\t' << fname;
  }
  file << '\n';

  // data
  for (size_t i=0; i<Npts; ++i) {
    file << lon[i] << '\t' << lat[i];
    for (size_t j=0; j<Nfld; ++j)
      file << '\t' << afvalues[j][i];
    file << '\n';
  }

  file.close();
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines


SimpleText* atlas__SimpleText__new()
{ return new SimpleText(); }


void atlas__SimpleText__delete (SimpleText* This)
{ delete This; }


grids::Unstructured* atlas__SimpleText__read (SimpleText* This, char* file_path)
{ return This->read(file_path); }


grids::Unstructured* atlas__read_simpletext (char* file_path)
{ return SimpleText::read(file_path); }


void atlas__write_simpletext_fieldset (char* file_path, FieldGroup* fieldset)
{ SimpleText::write(file_path, *fieldset); }


void atlas__write_simpletext_field (char* file_path, Field* field)
{ SimpleText::write(file_path, *field); }


// ------------------------------------------------------------------


} // namespace io
} // namespace atlas
