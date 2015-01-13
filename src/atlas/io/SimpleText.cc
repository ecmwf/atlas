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
#include <iterator>
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


template< typename DATA_TYPE >
void write_field_nodes(const Field& field, std::ostream& out)
{
#if 0
  ArrayView< DATA_TYPE > data(field);
#endif
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

  int nb_nodes;
  file >> nb_nodes;

#if 0
  if (!grid->has_function_space("nodes")) {
    std::vector< int > shape(2);
    shape[0] = nb_nodes;
    shape[1] = Field::UNDEF_VARS;
    grid->create_function_space("nodes","Lagrange_P0",shape)
        .metadata().set("type",static_cast< int >(Entity::NODES));
  }

  if (grid->function_space("nodes").shape(0) != nb_nodes) {
    throw Exception("existing nodes function space has incompatible number of nodes",Here());
  }

  FunctionSpace& nodes = grid->function_space("nodes");
  if (!nodes.has_field("coordinates"))
    nodes.create_field< double >("coordinates",3);

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
    const int npnt, const double*& lon, const double*& lat,
    const int nfld, const double** afields )
{
  const size_t
      Npnt (npnt>=0?         npnt : 0),
      Nfld (nfld>=0 && Npnt? nfld : 0);
  if (!Npnt)
    throw eckit::BadParameter("SimpleText::write: invalid number of points (npnt)");

  std::ofstream file(file_path.c_str());
  if (!file.is_open())
    throw eckit::CantOpenFile(file_path);

  // header
  file << Npnt << ' ' << Nfld << '\n';

  // data
  for (size_t i=0; i<Npnt; ++i) {
    file << lon[i] << '\t' << lon[i];
#if 0
    for (size_t j=0; j<Nfld; ++j)
      file << '\t' << afields[i][j];
#else
  std::copy(
    afields[i],
    afields[i] + sizeof(afields[i])/sizeof(double),
    std::ostream_iterator< double >(file,"\t") );
#endif
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
