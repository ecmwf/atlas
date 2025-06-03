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
#include <type_traits>
#include <vector>

#include "atlas/util/Object.h"

#include "atlas/library/config.h"

namespace eckit {
class Configuration;
}

namespace atlas {
class FieldSet;
class Field;
class Grid;
class Projection;
namespace util {
class Metadata;
class PartitionPolygon;
class PartitionPolygons;
}  // namespace util
namespace parallel {
class GatherScatter;
}  // namespace parallel

}  // namespace atlas

namespace atlas {
namespace functionspace {

#define FunctionspaceT_nonconst typename std::remove_const<FunctionSpaceT>::type
#define FunctionspaceT_const typename std::add_const<FunctionSpaceT>::type

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpaceImpl : public util::Object {
public:
    FunctionSpaceImpl();
    virtual ~FunctionSpaceImpl();
    virtual std::string type() const = 0;
    virtual operator bool() const { return true; }
    virtual size_t footprint() const = 0;

    virtual atlas::Field createField(const eckit::Configuration&) const = 0;

    virtual atlas::Field createField(const atlas::Field&, const eckit::Configuration&) const = 0;

    atlas::Field createField(const atlas::Field&) const;

    template <typename DATATYPE>
    atlas::Field createField(const eckit::Configuration&) const;

    template <typename DATATYPE>
    atlas::Field createField() const;

    const util::Metadata& metadata() const { return *metadata_; }
    util::Metadata& metadata() { return *metadata_; }

    template <typename FunctionSpaceT>
    FunctionspaceT_nonconst* cast();

    template <typename FunctionSpaceT>
    FunctionspaceT_const* cast() const;

    virtual std::string distribution() const = 0;

    virtual void haloExchange(const FieldSet&, bool /*on_device*/ = false) const;
    virtual void haloExchange(const Field&, bool /* on_device*/ = false) const;

    virtual void adjointHaloExchange(const FieldSet&, bool /*on_device*/ = false) const;
    virtual void adjointHaloExchange(const Field&, bool /* on_device*/ = false) const;

    virtual void gather(const FieldSet&, FieldSet&) const;
    virtual void gather(const Field&, Field&) const;

    virtual void scatter(const FieldSet&, FieldSet&) const;
    virtual void scatter(const Field&, Field&) const;

    virtual const parallel::GatherScatter& gather() const;
    virtual const parallel::GatherScatter& scatter() const;

    virtual idx_t size() const = 0;

    virtual idx_t part() const;

    virtual idx_t nb_parts() const;

    virtual const util::PartitionPolygon& polygon(idx_t halo = 0) const;

    virtual const atlas::Grid& grid() const;

    virtual atlas::Field lonlat() const;

    virtual atlas::Field ghost() const;

    virtual atlas::Field remote_index() const;

    virtual atlas::Field partition() const;

    virtual atlas::Field global_index() const;

    virtual const util::PartitionPolygons& polygons() const;

    virtual const Projection& projection() const;

    virtual std::string mpi_comm() const;

private:
    util::Metadata* metadata_;
};

template <typename FunctionSpaceT>
inline FunctionspaceT_nonconst* FunctionSpaceImpl::cast() {
    return dynamic_cast<FunctionspaceT_nonconst*>(this);
}

template <typename FunctionSpaceT>
inline FunctionspaceT_const* FunctionSpaceImpl::cast() const {
    return dynamic_cast<FunctionspaceT_const*>(this);
}

#undef FunctionspaceT_const
#undef FunctionspaceT_nonconst

//------------------------------------------------------------------------------------------------------

/// @brief Dummy Functionspace class that evaluates to false
class NoFunctionSpace : public FunctionSpaceImpl {
public:
    NoFunctionSpace(): FunctionSpaceImpl() {}
    virtual ~NoFunctionSpace() {}
    virtual std::string type() const { return "NoFunctionSpace"; }
    virtual operator bool() const { return false; }
    virtual size_t footprint() const { return sizeof(*this); }
    virtual std::string distribution() const { return std::string(); }

    virtual Field createField(const eckit::Configuration&) const;
    virtual Field createField(const Field&, const eckit::Configuration&) const;
    virtual idx_t size() const { return 0; }
};

//------------------------------------------------------------------------------------------------------

}  // namespace functionspace

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
