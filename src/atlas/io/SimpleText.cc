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
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/log/Log.h"
#include "atlas/Field.h"
#include "atlas/FieldGroup.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/util/ArrayView.h"

#include "atlas/io/SimpleText.h"


using namespace eckit;


namespace atlas {
namespace io {


namespace {


template< typename DATA_TYPE >
void write_field_nodes(const SimpleText& st, const Field& field, std::ostream& out)
{
  Log::info() << "writing field " << field.name() << "..." << std::endl;

#if 0
  ArrayView< DATA_TYPE > data(field);
#endif
}


} // end anonymous namespace


Mesh* SimpleText::read(const std::string& file_path, Mesh *mesh)
{
  if (mesh==NULL)
    mesh = new Mesh();

  std::ifstream file(file_path.c_str());
  if (!file.is_open())
    throw CantOpenFile(file_path);

  int nb_nodes;
  file >> nb_nodes;

  if (!mesh->has_function_space("nodes")) {
    std::vector< int > shape(2);
    shape[0] = nb_nodes;
    shape[1] = Field::UNDEF_VARS;
    mesh->create_function_space("nodes","Lagrange_P0",shape)
        .metadata().set("type",static_cast< int >(Entity::NODES));
  }

  if (mesh->function_space("nodes").shape(0) != nb_nodes) {
    throw Exception("existing nodes function space has incompatible number of nodes",Here());
  }

  FunctionSpace& nodes = mesh->function_space("nodes");
  if (!nodes.has_field("coordinates"))
    nodes.create_field< double >("coordinates",3);

  ArrayView< double, 2 > coords(nodes.field("coordinates"));
  for (size_t n=0; n<nb_nodes; ++n)
    file >> coords(n,XX) >> coords(n,YY);

  file.close();
  return mesh;
}


void SimpleText::write(const std::string& file_path, const Mesh &mesh) const
{
  LocalPathName path(file_path);
  Log::info() << "writing mesh to file " << path << std::endl;
  std::ofstream file(path.c_str());

  FunctionSpace& nodes = mesh.function_space( "nodes" );
  ArrayView< double, 2 > coords( nodes.field( "coordinates" ) );
  const int nb_nodes = nodes.shape(0);

  file << nb_nodes << '\n';
  for( size_t n = 0; n < nb_nodes; ++n )
    file << coords(n,XX) << ' ' << coords(n,YY) << '\n';

  file << std::flush;
  file.close();
}


void SimpleText::write(const std::string& file_path, const FieldGroup &fieldset) const
{
  LocalPathName path(file_path);
  Log::info() << "writing FieldGroup " << fieldset.name() << " to file " << path << std::endl;
  std::ofstream file(path.c_str());

  for (int i=0; i<fieldset.size(); ++i) {
    const Field& field = fieldset.field(i);
    if (field.function_space().metadata().has< int >("type") &&
        field.function_space().metadata().get< int >("type") == Entity::NODES) {
      field.data_type()=="int32"?  write_field_nodes< int    >(*this,field,file) :
      field.data_type()=="int64"?  write_field_nodes< long   >(*this,field,file) :
      field.data_type()=="real32"? write_field_nodes< float  >(*this,field,file) :
      field.data_type()=="real64"? write_field_nodes< double >(*this,field,file) :
                                   void();
    }
    else {
      throw Exception("function_space "+field.function_space().name()+" has no type.. ?");
    }
    file << std::flush;
  }

  file.close();
}


void SimpleText::write(const std::string& file_path, const Field &field) const
{
  LocalPathName path(file_path);
  Log::info() << "writing field " << field.name() << " to file " << path << std::endl;
  std::ofstream file(path.c_str());

  if (field.function_space().metadata().has< int >("type") &&
      field.function_space().metadata().get< int >("type") == Entity::NODES) {
    field.data_type()=="int32"?  write_field_nodes< int    >(*this,field,file) :
    field.data_type()=="int64"?  write_field_nodes< long   >(*this,field,file) :
    field.data_type()=="real32"? write_field_nodes< float  >(*this,field,file) :
    field.data_type()=="real64"? write_field_nodes< double >(*this,field,file) :
                                 void();
  }
  else {
    throw Exception("function_space "+field.function_space().name()+" has no type.. ?");
  }

  file << std::flush;
  file.close();
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines


SimpleText* atlas__SimpleText__new()
{ return new SimpleText(); }


void atlas__SimpleText__delete (SimpleText* This)
{ delete This; }


Mesh* atlas__SimpleText__read (SimpleText* This, char* file_path)
{ return This->read(file_path); }


void atlas__SimpleText__write (SimpleText* This, char* file_path, Mesh* mesh)
{ This->write(file_path, *mesh); }


Mesh* atlas__read_simpletext (char* file_path)
{ return SimpleText::read(file_path); }


void atlas__write_simpletext_mesh (char* file_path, Mesh* mesh) {
  SimpleText writer;
  writer.write(file_path, *mesh);
}


void atlas__write_simpletext_fieldset (char* file_path, FieldGroup* fieldset) {
  SimpleText writer;
  writer.write(file_path, *fieldset);
}


void atlas__write_simpletext_field (char* file_path, Field* field) {
  SimpleText writer;
  writer.write(file_path, *field);
}


// ------------------------------------------------------------------


} // namespace io
} // namespace atlas
