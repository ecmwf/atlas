/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Options.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/Earth.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

type::type(const std::string& _type) {
    set("type", _type);
}

halo::halo(size_t size) {
    set("halo", size);
}

datatype::datatype(DataType::kind_t kind) {
    set("datatype", kind);
}

datatype::datatype(const std::string& str) {
    set("datatype", DataType::str_to_kind(str));
}

datatype::datatype(DataType dtype) {
    set("datatype", dtype.kind());
}

name::name(const std::string& _name) {
    set("name", _name);
}

global::global(size_t _owner) {
    set("global", true);
    set("owner", _owner);
}

levels::levels(size_t _levels) {
    set("levels", _levels);
}

variables::variables(size_t _variables) {
    set("variables", _variables);
}

vector::vector(size_t _components) {
    set("variables", _components);
    set("type", "vector");
}

radius::radius(double _radius) {
    set("radius", _radius);
}

radius::radius(const std::string& key) {
    if (key == "Earth") {
        set("radius", util::Earth::radius());
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

pole_edges::pole_edges(bool _pole_edges) {
    set("pole_edges", _pole_edges);
}

alignment::alignment(int value) {
    set("alignment", value);
}

on_device::on_device() {
  set("execution_space", "device");
}

// ----------------------------------------------------------------------------

}  // namespace option
}  // namespace atlas
