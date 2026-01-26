/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/mesh/Halo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"

// ----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class FieldSet;
}

namespace atlas {
namespace parallel {
class HaloExchange;
class GatherScatter;
class Checksum;
}  // namespace parallel
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {

// ----------------------------------------------------------------------------

class CellColumns : public functionspace::FunctionSpaceImpl {
public:
    CellColumns(const Mesh&, const eckit::Configuration& = util::NoConfig());

    ~CellColumns() override;

    static std::string static_type() { return "CellColumns"; }
    std::string type() const override { return static_type(); }

    std::string distribution() const override;

    idx_t nb_cells() const;
    idx_t nb_cells_global() const;  // Only on MPI rank 0, will this be different from 0
    idx_t levels() const;
    std::vector<idx_t> nb_cells_global_foreach_rank() const;

    const Mesh& mesh() const { return mesh_; }
    Mesh& mesh() { return mesh_; }

    const mesh::HybridElements& cells() const { return cells_; }
    mesh::HybridElements& cells() { return cells_; }

    // -- Field creation methods

    Field createField(const eckit::Configuration&) const override;

    Field createField(const Field&, const eckit::Configuration&) const override;

    // -- Parallelisation aware methods

    const mesh::Halo& halo() const { return halo_; }

    void haloExchange(const FieldSet&, bool on_device = false) const override;
    void haloExchange(const Field&, bool on_device = false) const override;
    const parallel::HaloExchange& halo_exchange() const;

    void gather(const FieldSet&, FieldSet&) const override;
    void gather(const Field&, Field&) const override;
    const parallel::GatherScatter& gather() const override;

    void scatter(const FieldSet&, FieldSet&) const override;
    void scatter(const Field&, Field&) const override;
    const parallel::GatherScatter& scatter() const override;

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;
    const parallel::Checksum& checksum() const;

    idx_t size() const override { return nb_cells_; }

    const Grid& grid() const override;

    Field lonlat() const override;

    Field remote_index() const override;

    Field global_index() const override;

    Field ghost() const override;

    Field partition() const override;

    std::string mpi_comm() const override { return mesh_.mpi_comm(); }

private:  // methods
    idx_t config_size(const eckit::Configuration& config) const;
    array::DataType config_datatype(const eckit::Configuration&) const;
    std::string config_name(const eckit::Configuration&) const;
    idx_t config_levels(const eckit::Configuration&) const;
    array::ArrayShape config_shape(const eckit::Configuration&) const;
    void set_field_metadata(const eckit::Configuration&, Field&) const;
    size_t footprint() const override;

private:                           // data
    mutable Grid grid_;
    Mesh mesh_;                    // non-const because functionspace may modify mesh
    mesh::HybridElements& cells_;  // non-const because functionspace may modify mesh
    idx_t nb_levels_;
    mesh::Halo halo_;
    idx_t nb_cells_;
    mutable long nb_cells_global_{-1};
    mutable util::ObjectHandle<parallel::GatherScatter> gather_scatter_;  // without ghost
    mutable util::ObjectHandle<parallel::HaloExchange> halo_exchange_;
    mutable util::ObjectHandle<parallel::Checksum> checksum_;
};

// -------------------------------------------------------------------

extern "C" {

CellColumns* atlas__fs__CellColumns__new(Mesh::Implementation* mesh, const eckit::Configuration* config);
void atlas__fs__CellColumns__delete(CellColumns* This);
int atlas__fs__CellColumns__nb_cells(const CellColumns* This);
Mesh::Implementation* atlas__fs__CellColumns__mesh(CellColumns* This);
mesh::Cells* atlas__fs__CellColumns__cells(CellColumns* This);
field::FieldImpl* atlas__fs__CellColumns__create_field(const CellColumns* This, const eckit::Configuration* options);
field::FieldImpl* atlas__fs__CellColumns__create_field_template(const CellColumns* This,
                                                                const field::FieldImpl* field_template,
                                                                const eckit::Configuration* options);

void atlas__fs__CellColumns__halo_exchange_fieldset(const CellColumns* This, field::FieldSetImpl* fieldset);
void atlas__fs__CellColumns__halo_exchange_field(const CellColumns* This, field::FieldImpl* field);
const parallel::HaloExchange* atlas__fs__CellColumns__get_halo_exchange(const CellColumns* This);

void atlas__fs__CellColumns__gather_fieldset(const CellColumns* This, const field::FieldSetImpl* local,
                                             field::FieldSetImpl* global);
void atlas__fs__CellColumns__gather_field(const CellColumns* This, const field::FieldImpl* local,
                                          field::FieldImpl* global);
const parallel::GatherScatter* atlas__fs__CellColumns__get_gather(const CellColumns* This);

void atlas__fs__CellColumns__scatter_fieldset(const CellColumns* This, const field::FieldSetImpl* global,
                                              field::FieldSetImpl* local);
void atlas__fs__CellColumns__scatter_field(const CellColumns* This, const field::FieldImpl* global,
                                           field::FieldImpl* local);
const parallel::GatherScatter* atlas__fs__CellColumns__get_scatter(const CellColumns* This);

void atlas__fs__CellColumns__checksum_fieldset(const CellColumns* This, const field::FieldSetImpl* fieldset,
                                               char*& checksum, int& size, int& allocated);
void atlas__fs__CellColumns__checksum_field(const CellColumns* This, const field::FieldImpl* field, char*& checksum,
                                            int& size, int& allocated);
const parallel::Checksum* atlas__fs__CellColumns__get_checksum(const CellColumns* This);
}

}  // namespace detail

// -------------------------------------------------------------------

class CellColumns : public FunctionSpace {
public:
    using Implementation = detail::CellColumns;

public:
    CellColumns();
    CellColumns(const FunctionSpace&);
    CellColumns(const Mesh&, const eckit::Configuration&);
    CellColumns(const Mesh& mesh);

    static std::string type() { return detail::CellColumns::static_type(); }

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    idx_t nb_cells() const;
    idx_t nb_cells_global() const;  // Only on MPI rank 0, will this be different from 0
    idx_t levels() const;

    const Mesh& mesh() const;

    const mesh::HybridElements& cells() const;

    // -- Parallelisation aware methods
    const mesh::Halo& halo() const;

    const parallel::HaloExchange& halo_exchange() const;

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;
    const parallel::Checksum& checksum() const;

private:
    const detail::CellColumns* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas
