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

#include "atlas/functionspace/CellColumns.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace field {
class FieldSetImpl;
class FieldImpl;
}  // namespace field
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {

extern "C" {
const CellColumns* atlas__CellsFunctionSpace__new(Mesh::Implementation* mesh, const eckit::Configuration* config);
void atlas__CellsFunctionSpace__delete(CellColumns* This);
int atlas__CellsFunctionSpace__nb_cells(const CellColumns* This);
const Mesh::Implementation* atlas__CellsFunctionSpace__mesh(const CellColumns* This);
const mesh::HybridElements* atlas__CellsFunctionSpace__cells(const CellColumns* This);

void atlas__CellsFunctionSpace__halo_exchange_fieldset(const CellColumns* This, field::FieldSetImpl* fieldset);
void atlas__CellsFunctionSpace__halo_exchange_field(const CellColumns* This, field::FieldImpl* field);
const parallel::HaloExchange* atlas__CellsFunctionSpace__get_halo_exchange(const CellColumns* This);

void atlas__CellsFunctionSpace__gather_fieldset(const CellColumns* This, const field::FieldSetImpl* local,
                                                field::FieldSetImpl* global);
void atlas__CellsFunctionSpace__gather_field(const CellColumns* This, const field::FieldImpl* local,
                                             field::FieldImpl* global);
const parallel::GatherScatter* atlas__CellsFunctionSpace__get_gather(const CellColumns* This);

void atlas__CellsFunctionSpace__scatter_fieldset(const CellColumns* This, const field::FieldSetImpl* global,
                                                 field::FieldSetImpl* local);
void atlas__CellsFunctionSpace__scatter_field(const CellColumns* This, const field::FieldImpl* global,
                                              field::FieldImpl* local);
const parallel::GatherScatter* atlas__CellsFunctionSpace__get_scatter(const CellColumns* This);

void atlas__CellsFunctionSpace__checksum_fieldset(const CellColumns* This, const field::FieldSetImpl* fieldset,
                                                  char*& checksum, int& size, int& allocated);
void atlas__CellsFunctionSpace__checksum_field(const CellColumns* This, const field::FieldImpl* field, char*& checksum,
                                               int& size, int& allocated);
const parallel::Checksum* atlas__CellsFunctionSpace__get_checksum(const CellColumns* This);
}

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
