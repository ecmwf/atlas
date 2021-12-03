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

#include <vector>

//#include "atlas/util/Object.h"

#include "atlas/library/config.h"

namespace atlas {
class Vertical;
namespace field {
class FieldImpl;
}
}  // namespace atlas

namespace atlas {

extern "C" {
Vertical* atlas__Vertical__new(idx_t levels, const double z[]);
Vertical* atlas__Vertical__new_interval(idx_t levels, const double z[], const double interval[]);
void atlas__Vertical__delete(Vertical* This);
field::FieldImpl* atlas__Vertical__z(const Vertical* This);
int atlas__Vertical__size(const Vertical* This);
}

}  // namespace atlas
