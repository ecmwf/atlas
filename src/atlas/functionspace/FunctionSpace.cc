/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

namespace atlas {

FunctionSpace::FunctionSpace(): Handle(new functionspace::NoFunctionSpace()) {}


std::string FunctionSpace::type() const {
    return get()->type();
}

FunctionSpace::operator bool() const {
    return get()->operator bool();
}

size_t FunctionSpace::footprint() const {
    return get()->footprint();
}

const Grid& FunctionSpace::grid() const {
    return get()->grid();
}

Field FunctionSpace::createField(const eckit::Configuration& config) const {
    return get()->createField(config);
}

Field FunctionSpace::createField(const Field& other) const {
    return get()->createField(other);
}

Field FunctionSpace::createField(const Field& other, const eckit::Configuration& config) const {
    return get()->createField(other, config);
}

std::string FunctionSpace::distribution() const {
    return get()->distribution();
}

void FunctionSpace::haloExchange(const Field& field, bool on_device) const {
    get()->haloExchange(field, on_device);
}

void FunctionSpace::adjointHaloExchange(const Field& field, bool on_device) const {
    get()->adjointHaloExchange(field, on_device);
}

idx_t FunctionSpace::size() const {
    return get()->size();
}

idx_t FunctionSpace::part() const {
    return get()->part();
}

idx_t FunctionSpace::nb_parts() const {
    return get()->nb_parts();
}

Field FunctionSpace::lonlat() const {
    return get()->lonlat();
}

Field FunctionSpace::ghost() const {
    return get()->ghost();
}

Field FunctionSpace::global_index() const {
    return get()->global_index();
}

Field FunctionSpace::remote_index() const {
    return get()->remote_index();
}

Field FunctionSpace::partition() const {
    return get()->partition();
}

void FunctionSpace::haloExchange(const FieldSet& fields, bool on_device) const {
    get()->haloExchange(fields, on_device);
}

void FunctionSpace::adjointHaloExchange(const FieldSet& fields, bool on_device) const {
    get()->adjointHaloExchange(fields, on_device);
}

void FunctionSpace::gather(const FieldSet& local, FieldSet& global) const {
    get()->gather(local, global);
}

void FunctionSpace::gather(const Field& local, Field& global) const {
    get()->gather(local, global);
}

void FunctionSpace::scatter(const FieldSet& global, FieldSet& local) const {
    get()->scatter(global, local);
}

void FunctionSpace::scatter(const Field& global, Field& local) const {
    get()->scatter(global, local);
}

const parallel::GatherScatter& FunctionSpace::gather() const {
    return get()->gather();
}

const parallel::GatherScatter& FunctionSpace::scatter() const {
    return get()->scatter();
}

std::string FunctionSpace::mpi_comm() const {
    return get()->mpi_comm();
}

const util::PartitionPolygon& FunctionSpace::polygon(idx_t halo) const {
    return get()->polygon(halo);
}

const util::PartitionPolygons& FunctionSpace::polygons() const {
    return get()->polygons();
}

const Projection& FunctionSpace::projection() const {
    return get()->projection();
}

template <typename DATATYPE>
Field FunctionSpace::createField() const {
    return get()->createField<DATATYPE>();
}

template <typename DATATYPE>
Field FunctionSpace::createField(const eckit::Configuration& options) const {
    return get()->createField<DATATYPE>(options);
}

template Field FunctionSpace::createField<double>() const;
template Field FunctionSpace::createField<float>() const;
template Field FunctionSpace::createField<int>() const;
template Field FunctionSpace::createField<long>() const;

template Field FunctionSpace::createField<double>(const eckit::Configuration&) const;
template Field FunctionSpace::createField<float>(const eckit::Configuration&) const;
template Field FunctionSpace::createField<int>(const eckit::Configuration&) const;
template Field FunctionSpace::createField<long>(const eckit::Configuration&) const;


// ------------------------------------------------------------------

}  // namespace atlas
