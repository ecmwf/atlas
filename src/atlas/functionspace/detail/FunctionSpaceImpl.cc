/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "FunctionSpaceImpl.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/option/Options.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Metadata.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array/MakeView.h"

namespace atlas {
namespace functionspace {

// ------------------------------------------------------------------

FunctionSpaceImpl::FunctionSpaceImpl(): metadata_(new util::Metadata()) {}

FunctionSpaceImpl::~FunctionSpaceImpl() {
    delete metadata_;
}

atlas::Field FunctionSpaceImpl::createField(const atlas::Field& field) const {
    return createField(field, util::NoConfig());
}

void FunctionSpaceImpl::haloExchange(const FieldSet&, bool) const {
    throw_Exception("haloExchange() not implemented in derived class ["+type()+"]", Here());
}

void FunctionSpaceImpl::haloExchange(const Field&, bool) const {
    throw_Exception("haloExchange() not implemented in derived class ["+type()+"]", Here());
}

void FunctionSpaceImpl::adjointHaloExchange(const FieldSet&, bool) const {
    throw_Exception("adjointHaloExchange() not implemented in derived class ["+type()+"]", Here());
}

void FunctionSpaceImpl::adjointHaloExchange(const Field&, bool) const {
    throw_Exception("adjointHaloExchange() not implemented in derived class ["+type()+"]", Here());
}

Field NoFunctionSpace::createField(const eckit::Configuration&) const {
    throw_Exception("createField(const eckit::Configuration&) not implemented in derived class ["+type()+"]", Here());
}

Field NoFunctionSpace::createField(const Field&, const eckit::Configuration&) const {
    throw_Exception("createField(const Field&, const eckit::Configuration&) not implemented in derived class ["+type()+"]", Here());
}

const Grid& FunctionSpaceImpl::grid() const {
    throw_Exception("grid() not implemented in derived class ["+type()+"]", Here());
}

Field FunctionSpaceImpl::lonlat() const {
    throw_Exception("lonlat() not implemented in derived class ["+type()+"]", Here());
}

Field FunctionSpaceImpl::ghost() const {
    throw_Exception("ghost() not implemented in derived class ["+type()+"]", Here());
}

Field FunctionSpaceImpl::remote_index() const {
    throw_Exception("remote_index() not implemented in derived class ["+type()+"]", Here());
}

Field FunctionSpaceImpl::partition() const {
    throw_Exception("partition() not implemented in derived class ["+type()+"]", Here());
}

Field FunctionSpaceImpl::global_index() const {
    throw_Exception("global_index() not implemented in derived class ["+type()+"]", Here());
}

const util::PartitionPolygon& FunctionSpaceImpl::polygon(idx_t /*halo */) const {
    throw_Exception("polygon() not implemented in derived class", Here());
}

const util::PartitionPolygons& FunctionSpaceImpl::polygons() const {
    throw_Exception("polygons() not implemented in derived class", Here());
}

const Projection& FunctionSpaceImpl::projection() const {
    throw_Exception("projection() not implemented in derived class", Here());
}

template <typename DATATYPE>
Field FunctionSpaceImpl::createField(const eckit::Configuration& options) const {
    return createField(option::datatypeT<DATATYPE>() | options);
}

template <typename DATATYPE>
Field FunctionSpaceImpl::createField() const {
    return createField(option::datatypeT<DATATYPE>());
}

idx_t FunctionSpaceImpl::part() const {
    ATLAS_NOTIMPLEMENTED;
}

idx_t FunctionSpaceImpl::nb_parts() const {
    ATLAS_NOTIMPLEMENTED;
}

void FunctionSpaceImpl::gather(const FieldSet& local, FieldSet& global) const {
    ATLAS_NOTIMPLEMENTED;
}

void FunctionSpaceImpl::gather(const Field& local, Field& global) const {
    ATLAS_NOTIMPLEMENTED;
}

void FunctionSpaceImpl::scatter(const FieldSet& global, FieldSet& local) const {
    ATLAS_NOTIMPLEMENTED;
}

void FunctionSpaceImpl::scatter(const Field& global, Field& local) const {
    ATLAS_NOTIMPLEMENTED;
}

const parallel::GatherScatter& FunctionSpaceImpl::gather() const {
    ATLAS_NOTIMPLEMENTED;
}

const parallel::GatherScatter& FunctionSpaceImpl::scatter() const {
    ATLAS_NOTIMPLEMENTED;
}

std::string FunctionSpaceImpl::mpi_comm() const {
    return mpi::comm().name();
}

const HaloDescription& FunctionSpaceImpl::halo_description() const {
    if (not halo_description_) {
        halo_description_.reset(new HaloDescription(array::make_view<int,1>(ghost())));
    }
    return *halo_description_;
}

template Field FunctionSpaceImpl::createField<double>() const;
template Field FunctionSpaceImpl::createField<float>() const;
template Field FunctionSpaceImpl::createField<int>() const;
template Field FunctionSpaceImpl::createField<long>() const;


template Field FunctionSpaceImpl::createField<double>(const eckit::Configuration&) const;
template Field FunctionSpaceImpl::createField<float>(const eckit::Configuration&) const;
template Field FunctionSpaceImpl::createField<int>(const eckit::Configuration&) const;
template Field FunctionSpaceImpl::createField<long>(const eckit::Configuration&) const;


// ------------------------------------------------------------------

}  // namespace functionspace

// ------------------------------------------------------------------

}  // namespace atlas
