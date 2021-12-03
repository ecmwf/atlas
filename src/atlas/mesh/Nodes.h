/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   Nodes.h
/// @author Willem Deconinck
/// @date   August 2015

#pragma once

#include <map>
#include <string>

#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Topology.h"

namespace atlas {
namespace mesh {
template <class T>
class ConnectivityInterface;
class IrregularConnectivityImpl;
using IrregularConnectivity = ConnectivityInterface<IrregularConnectivityImpl>;
}  // namespace mesh
}  // namespace atlas

namespace atlas {
namespace mesh {

/**
 * \brief Nodes class that owns a collection of fields defined in nodes of the
 * mesh
 */
class Nodes : public util::Object {
public:
    using Connectivity = IrregularConnectivity;

    using Topology = util::Topology;

public:  // methods
         //-- Constructors
    /// @brief Construct "size" nodes
    Nodes();
    //  Nodes(idx_t size);

    //-- Accessors

    const Field& field(const std::string& name) const;
    Field& field(const std::string& name);
    bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }

    const Field& field(idx_t) const;
    Field& field(idx_t);
    idx_t nb_fields() const { return static_cast<idx_t>(fields_.size()); }

    const util::Metadata& metadata() const { return metadata_; }
    util::Metadata& metadata() { return metadata_; }

    const Field& global_index() const { return global_index_; }
    Field& global_index() { return global_index_; }

    const Field& remote_index() const { return remote_index_; }
    Field& remote_index() { return remote_index_; }

    const Field& partition() const { return partition_; }
    Field& partition() { return partition_; }

    const Field& xy() const { return xy_; }
    Field& xy() { return xy_; }

    const Field& lonlat() const { return lonlat_; }
    Field& lonlat() { return lonlat_; }

    const Field& ghost() const { return ghost_; }
    Field& ghost() { return ghost_; }

    const Field& flags() const { return flags_; }
    Field& flags() { return flags_; }

    const Field& halo() const { return halo_; }
    Field& halo() { return halo_; }

    /// @brief Node to Edge connectivity table
    const Connectivity& edge_connectivity() const;
    Connectivity& edge_connectivity();

    /// @brief Node to Cell connectivity table
    const Connectivity& cell_connectivity() const;
    Connectivity& cell_connectivity();

    const Connectivity& connectivity(const std::string& name) const;
    Connectivity& connectivity(const std::string& name);

    bool has_connectivity(std::string name) const { return connectivities_.count(name); }

    idx_t size() const { return size_; }

    // -- Modifiers

    Field add(const Field&);

    void resize(idx_t);

    void remove_field(const std::string& name);

    Connectivity& add(Connectivity*);

    /// @brief Return the memory footprint of the Nodes
    size_t footprint() const;

    void updateDevice() const;

    void updateHost() const;

    void syncHostDevice() const;

private:
    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Nodes& p) {
        p.print(s);
        return s;
    }

private:
    typedef std::map<std::string, Field> FieldMap;
    typedef std::map<std::string, util::ObjectHandle<Connectivity>> ConnectivityMap;

private:
    idx_t size_;
    FieldMap fields_;
    ConnectivityMap connectivities_;

    util::Metadata metadata_;

    // Cached shortcuts to specific fields in fields_
    Field global_index_;
    Field remote_index_;
    Field partition_;
    Field xy_;
    Field lonlat_;
    Field ghost_;
    Field flags_;
    Field halo_;

    // Cached shortcuts to specific connectivities in connectivities_
    Connectivity* edge_connectivity_;
    Connectivity* cell_connectivity_;
};

inline const Nodes::Connectivity& Nodes::edge_connectivity() const {
    return *edge_connectivity_;
}

inline Nodes::Connectivity& Nodes::edge_connectivity() {
    return *edge_connectivity_;
}

inline const Nodes::Connectivity& Nodes::cell_connectivity() const {
    return *cell_connectivity_;
}

inline Nodes::Connectivity& Nodes::cell_connectivity() {
    return *cell_connectivity_;
}

extern "C" {
Nodes* atlas__mesh__Nodes__create();
void atlas__mesh__Nodes__delete(Nodes* This);
idx_t atlas__mesh__Nodes__size(Nodes* This);
void atlas__mesh__Nodes__resize(Nodes* This, idx_t size);
idx_t atlas__mesh__Nodes__nb_fields(Nodes* This);
void atlas__mesh__Nodes__add_field(Nodes* This, field::FieldImpl* field);
void atlas__mesh__Nodes__remove_field(Nodes* This, char* name);
int atlas__mesh__Nodes__has_field(Nodes* This, char* name);
field::FieldImpl* atlas__mesh__Nodes__field_by_name(Nodes* This, char* name);
field::FieldImpl* atlas__mesh__Nodes__field_by_idx(Nodes* This, idx_t idx);
util::Metadata* atlas__mesh__Nodes__metadata(Nodes* This);
void atlas__mesh__Nodes__str(Nodes* This, char*& str, int& size);
IrregularConnectivity* atlas__mesh__Nodes__edge_connectivity(Nodes* This);
IrregularConnectivity* atlas__mesh__Nodes__cell_connectivity(Nodes* This);
IrregularConnectivity* atlas__mesh__Nodes__connectivity(Nodes* This, char* name);
void atlas__mesh__Nodes__add_connectivity(Nodes* This, IrregularConnectivity* connectivity);
field::FieldImpl* atlas__mesh__Nodes__xy(Nodes* This);
field::FieldImpl* atlas__mesh__Nodes__lonlat(Nodes* This);
field::FieldImpl* atlas__mesh__Nodes__global_index(Nodes* This);
field::FieldImpl* atlas__mesh__Nodes__remote_index(Nodes* This);
field::FieldImpl* atlas__mesh__Nodes__partition(Nodes* This);
field::FieldImpl* atlas__mesh__Nodes__ghost(Nodes* This);
}

//------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
