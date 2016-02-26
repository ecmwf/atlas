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
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "atlas/atlas_config.h"
#include "atlas/atlas_defines.h"
#include "atlas/grid/partitioners/Partitioner.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/util/runtime/Log.h"
#include "atlas/util/parallel/mpi/mpi.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/grid/partitioners/TransPartitioner.h"
#endif
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"


namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string,
           atlas::grid::partitioners::PartitionerFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string,
                         atlas::grid::partitioners::PartitionerFactory *>();
    }
}

namespace atlas {
namespace grid {
namespace partitioners {

Partitioner::Partitioner(const grid::Grid& grid): nb_partitions_(eckit::mpi::size()), grid_(grid)
{ }

Partitioner::Partitioner(const grid::Grid& grid, const size_t nb_partitions): nb_partitions_(nb_partitions), grid_(grid)
{ }

Partitioner::~Partitioner()
{ }

size_t Partitioner::nb_partitions() const
{
  return nb_partitions_;
}

GridDistribution* Partitioner::distribution() const
{
  return new GridDistribution(*this);
}


namespace {

template<typename T> void load_builder() { PartitionerBuilder<T>("tmp"); }

struct force_link {
    force_link()
    {
        load_builder< grid::partitioners::EqualRegionsPartitioner >();
#ifdef ATLAS_HAVE_TRANS
        load_builder< grid::partitioners::TransPartitioner >();
#endif
    }
};

}



PartitionerFactory::PartitionerFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


PartitionerFactory::~PartitionerFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void PartitionerFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, PartitionerFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

bool PartitionerFactory::has(const std::string& name)
{
  pthread_once(&once, init);

  eckit::AutoLock<eckit::Mutex> lock(local_mutex);

  static force_link static_linking;

  return ( m->find(name) != m->end() );
}



Partitioner* PartitionerFactory::build(const std::string& name, const grid::Grid& grid) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, PartitionerFactory *>::const_iterator j = m->find(name);

    Log::debug() << "Looking for PartitionerFactory [" << name << "]" << '\n';

    if (j == m->end()) {
        Log::error() << "No PartitionerFactory for [" << name << "]" << '\n';
        Log::error() << "PartitionerFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No PartitionerFactory called ") + name);
    }

    return (*j).second->make(grid);
}

Partitioner* PartitionerFactory::build(const std::string& name, const grid::Grid& grid, const size_t nb_partitions ) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, PartitionerFactory *>::const_iterator j = m->find(name);

    Log::debug() << "Looking for PartitionerFactory [" << name << "]" << '\n';

    if (j == m->end()) {
        Log::error() << "No PartitionerFactory for [" << name << "]" << '\n';
        Log::error() << "PartitionerFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No PartitionerFactory called ") + name);
    }

    return (*j).second->make(grid, nb_partitions);
}

} // namespace partitioners
} // namespace grid
} // namespace atlas
