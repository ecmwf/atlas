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

namespace eckit {
class Configuration;
}

namespace atlas {
namespace functionspace{
class FunctionSpace;
}
namespace field{
class Field;
class FieldSet;
}
}  // namespace atlas

namespace atlas {
namespace redistribution {
namespace detail {
class RedistributionImpl;

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {

RedistributionImpl* atlas__Redistribution__new__config(
    const functionspace::FunctionSpaceImpl* fspace1, const functionspace::FunctionSpaceImpl* fspace2,
    const eckit::Configuration* config);

void atlas__Redistribution__execute(
    const RedistributionImpl* This, const Field* field_1, Field* field_2);

const FunctionSpace* atlas__Redistribution__source(
    const RedistributionImpl* This);

const FunctionSpace* atlas__Redistribution__target(
    const RedistributionImpl* This);

}


}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
