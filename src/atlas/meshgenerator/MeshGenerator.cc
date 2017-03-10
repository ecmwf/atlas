/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/utils/MD5.h"

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/meshgenerator/DelaunayMeshGenerator.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"

using atlas::mesh::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, MeshGeneratorFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, MeshGeneratorFactory *>();
    }

    template<typename T> void load_builder() { MeshGeneratorBuilder<T>("tmp"); }

    struct force_link {
        force_link()
        {
            load_builder<meshgenerator::StructuredMeshGenerator>();
            load_builder<meshgenerator::DelaunayMeshGenerator>();
        }
    };

}

//----------------------------------------------------------------------------------------------------------------------

MeshGenerator* MeshGenerator::create(const std::string &key, const eckit::Parametrisation &params)
{
  return MeshGeneratorFactory::build(key,params);
}

MeshGenerator::MeshGenerator()
{
}

MeshGenerator::~MeshGenerator() {
}

Mesh* MeshGenerator::operator()( const grid::Grid& grid ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,*mesh);
  return mesh;
}

Mesh* MeshGenerator::operator()( const grid::Grid& grid, const grid::Distribution& distribution ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,distribution,*mesh);
  return mesh;
}

Mesh* MeshGenerator::generate( const grid::Grid& grid ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,*mesh);
  return mesh;
}

Mesh* MeshGenerator::generate( const grid::Grid& grid, const grid::Distribution& distribution ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,distribution,*mesh);
  return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

void MeshGenerator::generate_global_element_numbering( Mesh& mesh ) const
{
  size_t loc_nb_elems = mesh.cells().size();
  std::vector<size_t> elem_counts( parallel::mpi::comm().size() );
  std::vector<int> elem_displs( parallel::mpi::comm().size() );

  parallel::mpi::comm().allGather(loc_nb_elems, elem_counts.begin(), elem_counts.end());

  elem_displs.at(0) = 0;
  for(size_t jpart = 1; jpart < parallel::mpi::comm().size(); ++jpart)
  {
    elem_displs.at(jpart) = elem_displs.at(jpart-1) + elem_counts.at(jpart-1);
  }

  gidx_t gid = 1+elem_displs.at( parallel::mpi::comm().rank() );

  array::ArrayView<gidx_t,1> glb_idx( mesh.cells().global_index() );

  for( size_t jelem=0; jelem<mesh.cells().size(); ++jelem )
  {
    glb_idx(jelem) = gid++;
  }
}

//----------------------------------------------------------------------------------------------------------------------

MeshGeneratorFactory::MeshGeneratorFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


MeshGeneratorFactory::~MeshGeneratorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void MeshGeneratorFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, MeshGeneratorFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}


MeshGenerator *MeshGeneratorFactory::build(const std::string &name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, MeshGeneratorFactory *>::const_iterator j = m->find(name);

    Log::debug<ATLAS>() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No MeshGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "MeshGeneratorFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No MeshGeneratorFactory called ") + name);
    }

    return (*j).second->make();
}

MeshGenerator *MeshGeneratorFactory::build(const std::string& name, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, MeshGeneratorFactory *>::const_iterator j = m->find(name);

    Log::debug<ATLAS>() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No MeshGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "MeshGeneratorFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No MeshGeneratorFactory called ") + name);
    }

    return (*j).second->make(param);
}

//----------------------------------------------------------------------------------------------------------------------

extern "C" {

void atlas__MeshGenerator__delete(MeshGenerator* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

MeshGenerator* atlas__MeshGenerator__create_noconfig(const char* name)
{
  MeshGenerator* meshgenerator(0);
  ATLAS_ERROR_HANDLING (
    meshgenerator = MeshGenerator::create(std::string(name));
  );
  return meshgenerator;
}

MeshGenerator* atlas__MeshGenerator__create(const char* name, const eckit::Parametrisation* params)
{
  MeshGenerator* meshgenerator(0);
  ATLAS_ERROR_HANDLING (
    ASSERT(params);
    meshgenerator = MeshGenerator::create(std::string(name),*params);
  );
  return meshgenerator;
}

Mesh* atlas__MeshGenerator__generate__grid_griddist (const MeshGenerator* This, const grid::Grid::grid_t* grid, const grid::Distribution::impl_t* distribution )
{
  ATLAS_ERROR_HANDLING(
    return This->generate(grid::Grid(grid), grid::Distribution(distribution));
  );
  return 0;
}

Mesh* atlas__MeshGenerator__generate__grid (const MeshGenerator* This, const grid::Grid::grid_t* grid )
{
  ATLAS_ERROR_HANDLING(
    return This->generate(grid::Grid(grid));
  );
  return 0;
}

}

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerator
} // namespace atlas

