/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "CellColumnsInterface.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace functionspace {
namespace detail {

using atlas::FieldSet;
using atlas::field::FieldImpl;
using atlas::field::FieldSetImpl;

// ----------------------------------------------------------------------

extern "C" {
const CellColumns* atlas__CellsFunctionSpace__new(Mesh::Implementation* mesh, const eckit::Configuration* config) {
    ATLAS_ASSERT(mesh);
    Mesh m(mesh);
    return new CellColumns(m, *config);
}

void atlas__CellsFunctionSpace__delete(CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    delete (This);
}

int atlas__CellsFunctionSpace__nb_cells(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return This->nb_cells();
}

const Mesh::Implementation* atlas__CellsFunctionSpace__mesh(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return This->mesh().get();
}

const mesh::HybridElements* atlas__CellsFunctionSpace__cells(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->cells();
}

void atlas__CellsFunctionSpace__halo_exchange_fieldset(const CellColumns* This, field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    FieldSet f(fieldset);
    This->haloExchange(f);
}

void atlas__CellsFunctionSpace__halo_exchange_field(const CellColumns* This, field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    Field f(field);
    This->haloExchange(f);
}

const parallel::HaloExchange* atlas__CellsFunctionSpace__get_halo_exchange(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->halo_exchange();
}

void atlas__CellsFunctionSpace__gather_fieldset(const CellColumns* This, const field::FieldSetImpl* local,
                                                field::FieldSetImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet l(local);
    FieldSet g(global);
    This->gather(l, g);
}

void atlas__CellsFunctionSpace__gather_field(const CellColumns* This, const field::FieldImpl* local,
                                             field::FieldImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field l(local);
    Field g(global);
    This->gather(l, g);
}

const parallel::GatherScatter* atlas__CellsFunctionSpace__get_gather(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->gather();
}

const parallel::GatherScatter* atlas__CellsFunctionSpace__get_scatter(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->scatter();
}

void atlas__CellsFunctionSpace__scatter_fieldset(const CellColumns* This, const field::FieldSetImpl* global,
                                                 field::FieldSetImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet g(global);
    FieldSet l(local);
    This->scatter(g, l);
}

void atlas__CellsFunctionSpace__scatter_field(const CellColumns* This, const field::FieldImpl* global,
                                              field::FieldImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field g(global);
    Field l(local);
    This->scatter(g, l);
}

const parallel::Checksum* atlas__CellsFunctionSpace__get_checksum(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->checksum();
}

void atlas__CellsFunctionSpace__checksum_fieldset(const CellColumns* This, const field::FieldSetImpl* fieldset,
                                                  char*& checksum, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    std::string checksum_str(This->checksum(fieldset));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

void atlas__CellsFunctionSpace__checksum_field(const CellColumns* This, const field::FieldImpl* field, char*& checksum,
                                               int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::string checksum_str(This->checksum(field));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

}  // extern C

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
