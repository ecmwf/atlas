/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "eckit/types/FloatCompare.h"

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

using atlas::Grid;
using atlas::Projection;

namespace atlas {
namespace mesh {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

MeshImpl::MeshImpl(eckit::Stream&) {
    ATLAS_NOTIMPLEMENTED;
}

void MeshImpl::encode(eckit::Stream&) const {
    ATLAS_NOTIMPLEMENTED;
}

MeshImpl::MeshImpl(): nodes_(new mesh::Nodes()), dimensionality_(2) {
    metadata_.set("mpi_comm",mpi::comm().name());
    createElements();
}

MeshImpl::~MeshImpl() {
    while (mesh_observers_.size()) {
        MeshObserver* o = mesh_observers_.back();
        o->onMeshDestruction(*this);
        o->unregisterMesh(*this);  // will also delete observer from mesh
    }
}

void MeshImpl::print(std::ostream&) const {}

size_t MeshImpl::footprint() const {
    size_t size = sizeof(*this);

    size += metadata_.footprint();
    if (nodes_) {
        size += nodes_->footprint();
    }
    if (cells_) {
        size += cells_->footprint();
    }
    if (facets_) {
        size += facets_->footprint();
    }
    if (ridges_) {
        size += ridges_->footprint();
    }
    if (peaks_) {
        size += peaks_->footprint();
    }
    if (partition_graph_) {
        size += partition_graph_->footprint();
    }
    for (const auto& polygon : polygons_) {
        if (polygon) {
            size += polygon->footprint();
        }
    }

    return size;
}

void MeshImpl::createElements() {
    cells_.reset(new HybridElements());
    facets_.reset(new HybridElements());
    ridges_.reset(new HybridElements());
    peaks_.reset(new HybridElements());
    if (dimensionality_ == 2) {
        edges_ = facets_;
    }
    else if (dimensionality_ == 3) {
        edges_ = ridges_;
    }
    else {
        throw_Exception("Invalid Mesh dimensionality", Here());
    }

    ATLAS_ASSERT(edges_.owners() == 2);
}

bool MeshImpl::generated() const {
    return !(cells_->size() == 0 && facets_->size() == 0 && ridges_->size() == 0 && peaks_->size() == 0);
}

void MeshImpl::setProjection(const Projection& projection) {
    projection_ = projection;
}

void MeshImpl::setGrid(const Grid& grid) {
    grid_ = grid;
    if (not projection_) {
        projection_ = grid_.projection();
    }
}

idx_t MeshImpl::nb_partitions() const {
    idx_t n;
    if (not metadata().get("nb_parts", n)) {
        n = mpi::comm(mpi_comm()).size();
    }
    return n;
}

idx_t MeshImpl::partition() const {
    idx_t p;
    if (not metadata().get("part", p)) {
        p = mpi::comm(mpi_comm()).rank();
    }
    return p;
}

std::string MeshImpl::mpi_comm() const {
    return metadata().getString("mpi_comm");
}


void MeshImpl::updateDevice() const {
    if (nodes_) {
        nodes_->updateDevice();
    }
    if (cells_) {
        cells_->updateDevice();
    }
    if (facets_) {
        facets_->updateDevice();
    }
    if (ridges_) {
        ridges_->updateDevice();
    }
    if (peaks_) {
        peaks_->updateDevice();
    }
}

void MeshImpl::updateHost() const {
    if (nodes_) {
        nodes_->updateHost();
    }
    if (cells_) {
        cells_->updateHost();
    }
    if (facets_) {
        facets_->updateHost();
    }
    if (ridges_) {
        ridges_->updateHost();
    }
    if (peaks_) {
        peaks_->updateHost();
    }
}

void MeshImpl::syncHostDevice() const {
    if (nodes_) {
        nodes_->syncHostDevice();
    }
    if (cells_) {
        cells_->syncHostDevice();
    }
    if (facets_) {
        facets_->syncHostDevice();
    }
    if (ridges_) {
        ridges_->syncHostDevice();
    }
    if (peaks_) {
        peaks_->syncHostDevice();
    }
}

const PartitionGraph& MeshImpl::partitionGraph() const {
    if (not partition_graph_) {
        partition_graph_.reset(build_partition_graph(*this));
    }
    return *partition_graph_;
}

PartitionGraph::Neighbours MeshImpl::nearestNeighbourPartitions() const {
    return partitionGraph().nearestNeighbours(partition());
}

const PartitionPolygon& MeshImpl::polygon(idx_t halo) const {
    if (halo >= static_cast<idx_t>(polygons_.size())) {
        polygons_.resize(halo + 1);
    }
    if (not polygons_[halo]) {
        int mesh_halo = 0;
        metadata().get("halo", mesh_halo);
        if (halo > mesh_halo) {
            throw_Exception("Mesh does not contain a halo of size " + std::to_string(halo) + ".", Here());
        }

        polygons_[halo].reset(new PartitionPolygon(*this, halo));
    }
    return *polygons_[halo];
}

const util::PartitionPolygons& MeshImpl::polygons() const {
    if (all_polygons_.size() == 0) {
        polygon().allGather(all_polygons_);
    }
    return all_polygons_;
}

void MeshImpl::attachObserver(MeshObserver& observer) const {
    if (std::find(mesh_observers_.begin(), mesh_observers_.end(), &observer) == mesh_observers_.end()) {
        mesh_observers_.push_back(&observer);
    }
}

void MeshImpl::detachObserver(MeshObserver& observer) const {
    mesh_observers_.erase(std::remove(mesh_observers_.begin(), mesh_observers_.end(), &observer),
                          mesh_observers_.end());
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace mesh
}  // namespace atlas
