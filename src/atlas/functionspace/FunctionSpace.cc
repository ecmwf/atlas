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
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

namespace atlas {

FunctionSpace::FunctionSpace() : Handle( new functionspace::NoFunctionSpace() ) {}


std::string FunctionSpace::type() const {
    return get()->type();
}

FunctionSpace::operator bool() const {
    return get()->operator bool();
}

size_t FunctionSpace::footprint() const {
    return get()->footprint();
}

Field FunctionSpace::createField( const eckit::Configuration& config ) const {
    return get()->createField( config );
}

Field FunctionSpace::createField( const Field& other ) const {
    return get()->createField( other );
}

Field FunctionSpace::createField( const Field& other, const eckit::Configuration& config ) const {
    return get()->createField( other, config );
}

std::string FunctionSpace::distribution() const {
    return get()->distribution();
}

void FunctionSpace::haloExchange( const Field& field, bool on_device ) const {
    return get()->haloExchange( field, on_device );
}

void FunctionSpace::adjointHaloExchange( const Field& field, bool on_device ) const {
    return get()->adjointHaloExchange( field, on_device );
}


idx_t FunctionSpace::size() const {
    return get()->size();
}

idx_t FunctionSpace::nb_partitions() const {
    return get()->nb_partitions();
}

Field FunctionSpace::lonlat() const {
    return get()->lonlat();
}

Field FunctionSpace::ghost() const {
    return get()->ghost();
}

void FunctionSpace::haloExchange( const FieldSet& fields, bool on_device ) const {
    return get()->haloExchange( fields, on_device );
}

void FunctionSpace::adjointHaloExchange( const FieldSet& fields, bool on_device ) const {
    return get()->adjointHaloExchange( fields, on_device );
}

const util::PartitionPolygon& FunctionSpace::polygon( idx_t halo ) const {
    return get()->polygon( halo );
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
Field FunctionSpace::createField( const eckit::Configuration& options ) const {
    return get()->createField<DATATYPE>( options );
}

template Field FunctionSpace::createField<double>() const;
template Field FunctionSpace::createField<float>() const;
template Field FunctionSpace::createField<int>() const;
template Field FunctionSpace::createField<long>() const;

template Field FunctionSpace::createField<double>( const eckit::Configuration& ) const;
template Field FunctionSpace::createField<float>( const eckit::Configuration& ) const;
template Field FunctionSpace::createField<int>( const eckit::Configuration& ) const;
template Field FunctionSpace::createField<long>( const eckit::Configuration& ) const;


// ------------------------------------------------------------------

}  // namespace atlas
