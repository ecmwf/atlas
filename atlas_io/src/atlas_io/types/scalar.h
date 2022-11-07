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

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const int&, atlas::io::Metadata&);
size_t encode_metadata(const long&, atlas::io::Metadata&);
size_t encode_metadata(const long long&, atlas::io::Metadata&);
size_t encode_metadata(const unsigned long&, atlas::io::Metadata&);
size_t encode_metadata(const unsigned long long&, atlas::io::Metadata&);
size_t encode_metadata(const float&, atlas::io::Metadata&);
size_t encode_metadata(const double&, atlas::io::Metadata&);

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const int&, atlas::io::Data&);
void encode_data(const long&, atlas::io::Data&);
void encode_data(const long long&, atlas::io::Data&);
void encode_data(const unsigned long&, atlas::io::Data&);
void encode_data(const unsigned long long&, atlas::io::Data&);
void encode_data(const float&, atlas::io::Data&);
void encode_data(const double&, atlas::io::Data&);

//---------------------------------------------------------------------------------------------------------------------

void decode(const atlas::io::Metadata&, const atlas::io::Data&, int&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, long&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, long long&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, unsigned long&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, unsigned long long&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, float&);
void decode(const atlas::io::Metadata&, const atlas::io::Data&, double&);

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
