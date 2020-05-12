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

#include <string>

#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Configuration;
}

namespace atlas {
class Field;
class FieldSet;
class Projection;
namespace functionspace {
class FunctionSpaceImpl;
}
namespace util {
class PartitionPolygon;
class PartitionPolygons;
}  // namespace util
}  // namespace atlas

namespace atlas {

//------------------------------------------------------------------------------------------------------

class FunctionSpace : DOXYGEN_HIDE( public util::ObjectHandle<functionspace::FunctionSpaceImpl> ) {
public:
    using Handle::Handle;
    FunctionSpace();

    std::string type() const;
    operator bool() const;
    size_t footprint() const;
    std::string distribution() const;

    Field createField( const eckit::Configuration& ) const;

    Field createField( const Field& ) const;
    Field createField( const Field&, const eckit::Configuration& ) const;

    template <typename DATATYPE>
    Field createField( const eckit::Configuration& ) const;

    template <typename DATATYPE>
    Field createField() const;

    void haloExchange( const FieldSet&, bool on_device = false ) const;
    void haloExchange( const Field&, bool on_device = false ) const;

    void adjointHaloExchange( const FieldSet&, bool on_device = false ) const;
    void adjointHaloExchange( const Field&, bool on_device = false ) const;

    const util::PartitionPolygon& polygon( idx_t halo = 0 ) const;

    const util::PartitionPolygons& polygons() const;

    const Projection& projection() const;

    idx_t nb_partitions() const;

    idx_t size() const;

    Field lonlat() const;

    Field ghost() const;
};

//------------------------------------------------------------------------------------------------------

extern template Field FunctionSpace::createField<float>() const;
extern template Field FunctionSpace::createField<double>() const;
extern template Field FunctionSpace::createField<int>() const;
extern template Field FunctionSpace::createField<long>() const;
extern template Field FunctionSpace::createField<float>( const eckit::Configuration& ) const;
extern template Field FunctionSpace::createField<double>( const eckit::Configuration& ) const;
extern template Field FunctionSpace::createField<int>( const eckit::Configuration& ) const;
extern template Field FunctionSpace::createField<long>( const eckit::Configuration& ) const;

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
