/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/detail/Cache.h"

namespace atlas {
namespace functionspace {


// Helper functions to get fields.
namespace {
template <typename BaseFunctionSpace>
Field getTij(const Mesh& mesh);

template <typename BaseFunctionSpace>
Field getGhost(const Mesh& mesh);

template <>
Field getTij<NodeColumns>(const Mesh& mesh) {
    return mesh.nodes().field("tij");
}

template <>
Field getTij<CellColumns>(const Mesh& mesh) {
    return mesh.cells().field("tij");
}

template <>
Field getGhost<NodeColumns>(const Mesh& mesh) {
    return mesh.nodes().ghost();
}

template <>
Field getGhost<CellColumns>(const Mesh& mesh) {
    // No ghost field for CellColumns. Halo field is next best thing.
    return mesh.cells().halo();
}
}  // namespace

namespace {
template <typename BaseFunctionSpace>
class CubedSphereStructureCache : public util::Cache<std::string, detail::CubedSphereStructure>,
                                  public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, detail::CubedSphereStructure>;
    CubedSphereStructureCache(): Base("CubedSphereStructureCache<" + BaseFunctionSpace::type() + ">") {}

public:
    static CubedSphereStructureCache& instance() {
        static CubedSphereStructureCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const BaseFunctionSpace* functionspace) {
        ATLAS_ASSERT(functionspace);
        auto& mesh = functionspace->mesh();
        ATLAS_ASSERT(mesh);
        auto& mesh_impl = *mesh.get();
        registerMesh(mesh_impl);
        creator_type creator = std::bind(&CubedSphereStructureCache::create, mesh, functionspace->size());
        util::ObjectHandle<value_type> value = Base::get_or_create(key(mesh_impl), creator);
        return value;
    }
    void onMeshDestruction(mesh::detail::MeshImpl& mesh) override { remove(key(mesh)); }

private:
    static Base::key_type key(const mesh::detail::MeshImpl& mesh) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "]";
        return key.str();
    }

    static value_type* create(const Mesh& mesh, idx_t size) {
        value_type* value = new value_type(getTij<BaseFunctionSpace>(mesh), getGhost<BaseFunctionSpace>(mesh), size);
        return value;
    }
};


}  // namespace

// All constructors pass arguments through to BaseFunctionSpace, then construct
// CubedSphereStructure.
template <typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns():
    BaseFunctionSpace(), cubedSphereColumnsHandle_(new detail::CubedSphereStructure()) {}

template <typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(const FunctionSpace& functionspace):
    BaseFunctionSpace([&]() {
        bool compatible = dynamic_cast<const typename BaseFunctionSpace::Implementation*>(functionspace.get());
        if (not compatible) {
            ATLAS_THROW_EXCEPTION("FunctionSpace " << functionspace.type() << " can not be interpreted as a "
                                                   << BaseFunctionSpace::type());
        }
        return functionspace;
    }()),
    cubedSphereColumnsHandle_(CubedSphereStructureCache<BaseFunctionSpace>::instance().get_or_create(this)) {}

template <typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(const Mesh& mesh, const eckit::Configuration& configuration):
    BaseFunctionSpace(mesh, configuration),
    cubedSphereColumnsHandle_(CubedSphereStructureCache<BaseFunctionSpace>::instance().get_or_create(this)) {}

template <typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(const Mesh& mesh):
    BaseFunctionSpace(mesh),
    cubedSphereColumnsHandle_(CubedSphereStructureCache<BaseFunctionSpace>::instance().get_or_create(this)) {}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::invalid_index() const {
    return cubedSphereColumnsHandle_.get()->invalid_index();
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::sizeOwned() const {
    return cubedSphereColumnsHandle_.get()->sizeOwned();
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::i_begin(idx_t t) const {
    return cubedSphereColumnsHandle_.get()->i_begin(t);
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::i_end(idx_t t) const {
    return cubedSphereColumnsHandle_.get()->i_end(t);
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::j_begin(idx_t t) const {
    return cubedSphereColumnsHandle_.get()->j_begin(t);
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::j_end(idx_t t) const {
    return cubedSphereColumnsHandle_.get()->j_end(t);
}

template <typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::index(idx_t t, idx_t i, idx_t j) const {
    return cubedSphereColumnsHandle_.get()->index(t, i, j);
}

template <typename BaseFunctionSpace>
bool CubedSphereColumns<BaseFunctionSpace>::is_valid_index(idx_t t, idx_t i, idx_t j) const {
    return cubedSphereColumnsHandle_.get()->is_valid_index(t, i, j);
}

template <typename BaseFunctionSpace>
Field CubedSphereColumns<BaseFunctionSpace>::tij() const {
    return cubedSphereColumnsHandle_.get()->tij();
}

// Explicit instantiation of template classes.
template class CubedSphereColumns<CellColumns>;
template class CubedSphereColumns<NodeColumns>;

}  // namespace functionspace
}  // namespace atlas
