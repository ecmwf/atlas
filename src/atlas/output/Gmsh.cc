/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/io/Gmsh.h"

using atlas::field::Field;
using atlas::field::FieldSet;
using atlas::mesh::Mesh;
using atlas::functionspace::FunctionSpace;
using eckit::Parametrisation;

namespace atlas {
namespace output {

void Gmsh::defaults()
{
  config_.binary = false;
  config_.nodes = "lonlat";
  config_.gather = false;
  config_.ghost = false;
  config_.elements = true;
  config_.edges = false;
  config_.radius = 1.0;
  config_.levels.clear();
  config_.file = "output.msh";
  config_.info = false;
  config_.openmode = "w";
}

namespace {
void merge(Gmsh::Configuration& present, const eckit::Parametrisation& update)
{
  update.get("binary",present.binary);
  update.get("nodes",present.nodes);
  update.get("gather",present.gather);
  update.get("ghost",present.ghost);
  update.get("elements",present.elements);
  update.get("edges",present.edges);
  update.get("radius",present.radius);
  update.get("levels",present.levels);
  update.get("file",present.file);
  update.get("info",present.info);
  update.get("openmode",present.openmode);
}
util::io::Gmsh writer(const Gmsh::Configuration& c)
{
  util::io::Gmsh gmsh;
  gmsh.options.set("ascii", not c.binary);
  gmsh.options.set("nodes",c.nodes);
  gmsh.options.set("gather",c.gather);
  gmsh.options.set("ghost",c.ghost);
  gmsh.options.set("elements",c.elements);
  gmsh.options.set("edges",c.edges);
  gmsh.options.set("radius",c.radius);
  gmsh.options.set("levels",c.levels);
  gmsh.options.set("info",c.info);
  return gmsh;
}
std::ios_base::openmode openmode(const Gmsh::Configuration& c)
{
  std::ios_base::openmode omode;
  if     ( std::string(c.openmode)=="w" )  omode = std::ios_base::out;
  else if( std::string(c.openmode)=="a" )  omode = std::ios_base::app;
  return omode;
}

}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(Stream& stream)
{
  defaults();
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(Stream& stream,const eckit::Parametrisation& config)
{
  defaults();
  merge(config_,config);
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(const PathName& file, const std::string& mode)
{
  defaults();
  config_.file = file.asString();
  config_.openmode = std::string(mode);
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(const PathName& file, const std::string& mode, const eckit::Parametrisation& config)
{
  defaults();
  merge(config_,config);
  config_.file = file.asString();
  config_.openmode = std::string(mode);
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(const PathName& file)
{
  defaults();
  config_.file = file.asString();
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh(const PathName& file, const eckit::Parametrisation& config)
{
  defaults();
  merge(config_,config);
  config_.file = file.asString();
}

// -----------------------------------------------------------------------------

Gmsh::~Gmsh()
{
}

// -----------------------------------------------------------------------------

void Gmsh::output(
        const mesh::Mesh& mesh,
        const eckit::Parametrisation& config) const
{
  Gmsh::Configuration c = config_;
  merge(c,config);
  writer(c).write(mesh,c.file);
}

// -----------------------------------------------------------------------------

void Gmsh::output(
    const field::Field& field,
    const eckit::Parametrisation& config ) const
{
  Gmsh::Configuration c = config_;
  merge(c,config);
  writer(c).write(field,c.file,openmode(c));
}

// -----------------------------------------------------------------------------

void Gmsh::output(
    const field::FieldSet& fields,
    const eckit::Parametrisation& config) const
{
  Gmsh::Configuration c = config_;
  merge(c,config);
  writer(c).write(fields,fields.field(0).functionspace(),c.file,openmode(c));
}

// -----------------------------------------------------------------------------

void Gmsh::output(
    const field::Field& field,
    const functionspace::FunctionSpace& functionspace,
    const eckit::Parametrisation& config) const
{
  Gmsh::Configuration c = config_;
  merge(c,config);
  writer(c).write(field,functionspace,c.file,openmode(c));
}

// -----------------------------------------------------------------------------

void Gmsh::output(
    const field::FieldSet& fields,
    const functionspace::FunctionSpace& functionspace,
    const eckit::Parametrisation& config) const
{
  Gmsh::Configuration c = config_;
  merge(c,config);
  writer(c).write(fields,functionspace,c.file,openmode(c));
}

// -----------------------------------------------------------------------------

namespace {
static OutputBuilder< Gmsh > __gmsh("gmsh");
}

} // namespace output
} // namespace atlas

