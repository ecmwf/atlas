/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "PartitionerFromPrePartitionedMesh.h"

#include <map>
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Once.h"
#include "eckit/thread/Mutex.h"
#include "eckit/exception/Exceptions.h"


namespace atlas {
namespace grid {
namespace partitioners {


#if 0
namespace {


typedef std::map<std::string, PartitionerFromPrePartitionedMeshFactory*> PartitionerFromPrePartitionedMeshFactoryMap_t;
static PartitionerFromPrePartitionedMeshFactoryMap_t *m = 0;
static eckit::Mutex *local_mutex = 0;
static pthread_once_t once = PTHREAD_ONCE_INIT;


static void init() {
    local_mutex = new eckit::Mutex();
    m = new PartitionerFromPrePartitionedMeshFactoryMap_t();
}


}  // (anonymous namespace)


PartitionerFromPrePartitionedMeshFactory::PartitionerFromPrePartitionedMeshFactory(const std::string& name):
    name_(name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    if (m->find(name) != m->end()) {
        throw eckit::SeriousBug("PartitionerFromPrePartitionedMeshFactory duplicate '" + name + "'");
    }

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


PartitionerFromPrePartitionedMeshFactory::~PartitionerFromPrePartitionedMeshFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


PartitionerFromPrePartitionedMesh* PartitionerFromPrePartitionedMeshFactory::build(const std::string& name, const PartitionerFromPrePartitionedMesh::Config& config) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    PartitionerFromPrePartitionedMeshFactoryMap_t::const_iterator j = m->find(name);
    if (j == m->end()) {
        eckit::Log::error() << "PartitionerFromPrePartitionedMeshFactory '" << name << "' not found." << std::endl;
        eckit::Log::error() << "PartitionerFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j) {
            eckit::Log::error() << '\t' << (*j).first << std::endl;
        }
        throw eckit::SeriousBug("PartitionerFromPrePartitionedMeshFactory '" + name + "' not found.");
    }

    return (*j).second->make(config);
}
#endif


void PartitionerFromPrePartitionedMesh::dump() {
#if 0
    std::vector<double> x,y;
    x.reserve(node_partition.size());
    y.reserve(node_partition.size());
    for (size_t i=0; i<node_partition.size(); ++i) {
        if (node_partition[i] == mpi_rank) {
            x.push_back(lonlat_tgt_pts[i].lon());
            y.push_back(lonlat_tgt_pts[i].lat());
        }
    }
    size_t count = x.size();
    size_t count_all = x.size();
    comm.allReduceInPlace(count_all, eckit::mpi::sum());

    for (int r = 0; r < comm.size(); ++r) {
        if (mpi_rank == r) {
            std::ofstream f("partition.py", mpi_rank == 0? std::ios::trunc : std::ios::app);

            if (mpi_rank == 0) {
                f << "\n" "import matplotlib.pyplot as plt"
                     "\n" ""
                     "\n" "from itertools import cycle"
                     "\n" "cycol = cycle('bgrcmy').next"
                     "\n" ""
                     "\n" "fig = plt.figure()"
                     "\n" "ax = fig.add_subplot(111)"
                     "\n" "";
            }
            f << "\n" "count_" << r << " = " << count <<
                 "\n" "count_all_" << r << " = " << count_all <<
                 "\n" ""
                 "\n" "x_" << r << " = ["; for (std::vector<double>::const_iterator ix=x.begin(); ix!=x.end(); ++ix) { f << *ix << ", "; } f << "]"
                 "\n" "y_" << r << " = ["; for (std::vector<double>::const_iterator iy=y.begin(); iy!=y.end(); ++iy) { f << *iy << ", "; } f << "]"
                 "\n"
                 "\n" "c = cycol()"
                 "\n" "ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')"
                 "\n" "";
            if (mpi_rank == int(comm.size()) - 1) {
                f << "\n" "ax.set_xlim(  0-5, 360+5)"
                     "\n" "ax.set_ylim(-90-5,  90+5)"
                     "\n" "plt.show()";
            }
        }
        comm.barrier();
    }
//    exit(0);
#endif
}


}  // partitioners
}  // grid
}  // atlas

