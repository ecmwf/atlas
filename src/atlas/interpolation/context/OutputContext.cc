/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/context/OutputContext.h"

#include <map>
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/Once.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/grid/partitioners/PartitionerFromPrePartitionedMesh.h"
#include "atlas/grid/partitioners/Partitioner.h"
#include "atlas/grid/Structured.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"


namespace atlas {
namespace interpolation {
namespace context {


namespace {


typedef std::map<std::string, OutputContextFactory*> OutputContextFactoryMap_t;
static OutputContextFactoryMap_t *m = 0;
static eckit::Mutex *local_mutex = 0;
static pthread_once_t once = PTHREAD_ONCE_INIT;


static void init() {
    local_mutex = new eckit::Mutex();
    m = new OutputContextFactoryMap_t();
}


}  // (anonymous namespace)


OutputContext::OutputContext(
        const std::string& gridname,
        const std::string& partitioner,
        const std::string& meshGenerator,
        const mesh::Mesh::Ptr prePartitionedMesh,
        const grid::Domain& prePartitionedDomain,
        bool meshGeneratorTriangulate,
        double meshGeneratorAngle ) :
    Context(gridname, partitioner, meshGenerator, meshGeneratorTriangulate, meshGeneratorAngle) {
    using grid::partitioners::Partitioner;
    using mesh::generators::MeshGenerator;

    grid_.reset(grid::Structured::create(optionGridname_));
    ASSERT(grid_);

    grid::partitioners::Partitioner::Ptr partitioner_(grid::partitioners::PartitionerFactory::build(optionPartitioner_, *grid_));
    ASSERT(partitioner_);

    try {
        grid::partitioners::PartitionerFromPrePartitionedMesh& partner = dynamic_cast< grid::partitioners::PartitionerFromPrePartitionedMesh& >(*partitioner_);
        partner.setup(prePartitionedMesh, prePartitionedDomain);
        distribution_.reset(partner.distribution());
    } catch (std::bad_cast&) {
        throw eckit::UserError("Partitioner has to be a PartitionerFromPrePartitionedMesh");
    }
    ASSERT(distribution_);

    MeshGenerator::Ptr meshgen(mesh::generators::MeshGeneratorFactory::build(optionMeshGenerator_, meshGeneratorParams_));
    mesh_.reset(meshgen->generate(*grid_, *distribution_));
}


OutputContextFactory::OutputContextFactory(const std::string& name):
    name_(name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    if (m->find(name) != m->end()) {
        throw eckit::SeriousBug("InterpolationFactory duplicate '" + name + "'");
    }

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


OutputContextFactory::~OutputContextFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


OutputContext* OutputContextFactory::build(const std::string& key, const std::string& name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    OutputContextFactoryMap_t::const_iterator j = m->find(key);
    if (j == m->end()) {
        eckit::Log::error() << "InterpolationFactory '" << key << "' not found." << std::endl;
        eckit::Log::error() << "InterpolationFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j) {
            eckit::Log::error() << '\t' << (*j).first << std::endl;
        }
        throw eckit::SeriousBug("InterpolationFactory '" + key + "' not found.");
    }

    return (*j).second->make(name);
}


}  // context
}  // interpolation
}  // atlas
