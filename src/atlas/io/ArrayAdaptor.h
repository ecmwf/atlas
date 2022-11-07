/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

namespace atlas {
namespace io {
class Metadata;
class Data;
class ArrayReference;
}  // namespace io
}  // namespace atlas

namespace atlas {
namespace array {

class Array;

//---------------------------------------------------------------------------------------------------------------------

void interprete(const Array&, atlas::io::ArrayReference&);

//---------------------------------------------------------------------------------------------------------------------

void decode(const atlas::io::Metadata&, const atlas::io::Data&, Array&);

//---------------------------------------------------------------------------------------------------------------------
}  // namespace array
}  // namespace atlas
