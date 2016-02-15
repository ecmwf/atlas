/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/meshgen/MeshGenerator.h"

#include <string>
#include <map>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/Mesh.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/meshgen/Delaunay.h"

namespace atlas {
namespace meshgen {

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
            load_builder<ReducedGridMeshGenerator>();
            load_builder<Delaunay>();
        }
    };

}

MeshGenerator* MeshGenerator::create(const std::string &key, const eckit::Parametrisation &params)
{
  return MeshGeneratorFactory::build(key,params);
}

MeshGenerator::MeshGenerator()
{
}

MeshGenerator::~MeshGenerator() {
}

Mesh* MeshGenerator::operator()( const Grid& grid ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,*mesh);
  return mesh;
}

Mesh* MeshGenerator::operator()( const Grid& grid, const GridDistribution& distribution ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,distribution,*mesh);
  return mesh;
}

Mesh* MeshGenerator::generate( const Grid& grid ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,*mesh);
  return mesh;
}

Mesh* MeshGenerator::generate( const Grid& grid, const GridDistribution& distribution ) const
{
  Mesh* mesh = new Mesh;
  generate(grid,distribution,*mesh);
  return mesh;
}


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

    Log::info() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

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

    Log::info() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No MeshGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "MeshGeneratorFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No MeshGeneratorFactory called ") + name);
    }

    return (*j).second->make(param);
}

extern "C" {

void atlas__MeshGenerator__delete(MeshGenerator* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
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

Mesh* atlas__MeshGenerator__generate__grid_griddist (const MeshGenerator* This, const Grid* grid, const GridDistribution* distribution )
{
  ATLAS_ERROR_HANDLING(
    return This->generate(*grid,*distribution);
  );
  return 0;
}

Mesh* atlas__MeshGenerator__generate__grid (const MeshGenerator* This, const Grid* grid )
{
  ATLAS_ERROR_HANDLING(
    return This->generate(*grid);
  );
  return 0;
}


}

} // namespace meshgen
} // namespace atlas

