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
#include "atlas/option/Options.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Metadata.h"

namespace atlas {
namespace functionspace {

// ------------------------------------------------------------------

FunctionSpaceImpl::FunctionSpaceImpl() : metadata_( new util::Metadata() ) {}

FunctionSpaceImpl::~FunctionSpaceImpl() {
    delete metadata_;
}

atlas::Field FunctionSpaceImpl::createField( const atlas::Field& field ) const {
    return createField( field, util::NoConfig() );
}

void FunctionSpaceImpl::haloExchange( const FieldSet&, bool ) const {
    ATLAS_NOTIMPLEMENTED;
}

void FunctionSpaceImpl::haloExchange( const Field&, bool ) const {
    ATLAS_NOTIMPLEMENTED;
}

Field NoFunctionSpace::createField( const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
}

Field NoFunctionSpace::createField( const Field&, const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
}

const util::PartitionPolygon& FunctionSpaceImpl::polygon( idx_t halo ) const {
    throw_Exception( "polygon() not implemented in derived class", Here() );
}

const std::vector<util::PartitionPolygon*>& FunctionSpaceImpl::polygons() const {
    throw_Exception( "polygons() not implemented in derived class", Here() );
}

template <typename DATATYPE>
Field FunctionSpaceImpl::createField( const eckit::Configuration& options ) const {
    return createField( option::datatypeT<DATATYPE>() | options );
}

template <typename DATATYPE>
Field FunctionSpaceImpl::createField() const {
    return createField( option::datatypeT<DATATYPE>() );
}

idx_t FunctionSpaceImpl::nb_partitions() const {
    ATLAS_NOTIMPLEMENTED;
}


template Field FunctionSpaceImpl::createField<double>() const;
template Field FunctionSpaceImpl::createField<float>() const;
template Field FunctionSpaceImpl::createField<int>() const;
template Field FunctionSpaceImpl::createField<long>() const;


template Field FunctionSpaceImpl::createField<double>( const eckit::Configuration& ) const;
template Field FunctionSpaceImpl::createField<float>( const eckit::Configuration& ) const;
template Field FunctionSpaceImpl::createField<int>( const eckit::Configuration& ) const;
template Field FunctionSpaceImpl::createField<long>( const eckit::Configuration& ) const;


// ------------------------------------------------------------------

}  // namespace functionspace

// ------------------------------------------------------------------

}  // namespace atlas
