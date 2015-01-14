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
#include <algorithm>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/log/Log.h"
#include "atlas/Field.h"
#include "atlas/FieldGroup.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Grid.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/util/ArrayView.h"
#include "atlas/io/SimpleText.h"


//using namespace eckit;
using eckit::LocalPathName;
using eckit::Log;


namespace atlas {
namespace io {


namespace {


// versioning of input/output files to (maybe) help in the future
static const std::string SIMPLETEXT_VERSION("ST/1");


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


grids::Unstructured* SimpleText::read(const std::string& file_path)
{
  std::vector< Grid::Point >* pts;

  grids::Unstructured* grid = new grids::Unstructured(pts);
  ASSERT(grid);

  std::ifstream file(file_path.c_str());
  if (!file.is_open())
    throw eckit::CantOpenFile(file_path);

  size_t
      nb_nodes,
      nb_fields;
  file >> nb_nodes >> nb_fields;

  std::vector< int > shape(2);
  shape[0] = nb_nodes;
  shape[1] = Field::UNDEF_VARS;
  grid->mesh().create_function_space("nodes","Lagrange_P0",shape).metadata().set("type",static_cast< int >(Entity::NODES));

  FunctionSpace& nodes = grid->mesh().function_space("nodes");
  nodes.create_field< double >("coordinates",3);

#if 0
  ArrayView< double, 2 > coords(nodes.field("coordinates"));
  for (size_t n=0; n<nb_nodes; ++n)
    file >> coords(n,XX) >> coords(n,YY);
#endif

  file.close();
  return grid;
}


void SimpleText::write(const std::string &file_path, const grids::Unstructured &grid)
{
  NOTIMP;
}


void SimpleText::write(const std::string& file_path, const FieldGroup &fieldset)
{
  LocalPathName path(file_path);
  Log::info() << "writing FieldGroup " << fieldset.name() << " to file " << path << std::endl;
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
      throw eckit::Exception("function_space "+field.function_space().name()+" has no type.. ?");
    }
    file << std::flush;
  }

  file.close();
}


void SimpleText::write(const std::string& file_path, const Field &field)
{
  LocalPathName path(file_path);
  Log::info() << "writing field " << field.name() << " to file " << path << std::endl;
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
    throw eckit::Exception("function_space "+field.function_space().name()+" has no type.. ?");
  }

  file << std::flush;
  file.close();
}


void SimpleText::write(
    const std::string& file_path,
    const int& nb_nodes, const double*& lon, const double*& lat,
    const int& nb_fields, const double** afvalues, const char** afnames )
{
  const size_t
      Nnod (nb_nodes>0?  nb_nodes  : 0),
      Nfld (nb_fields>0? nb_fields : 0);
  if (!Nnod)
    throw eckit::BadParameter("SimpleText::write: invalid number of points (n)");

  std::ofstream file(file_path.c_str());
  if (!file.is_open())
    throw eckit::CantOpenFile(file_path);

  // header
  file << SIMPLETEXT_VERSION << ' ' << Nnod << ' ' << Nfld << "\n"
          "lon\tlat";
  for (size_t j=0; j<Nfld; ++j) {
    std::string fname(afnames[j]);
    sanitize_field_name(fname);
    file << '\t' << fname;
  }
  file << '\n';

  // data
  for (size_t i=0; i<Nnod; ++i) {
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
