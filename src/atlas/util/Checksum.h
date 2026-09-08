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

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>

namespace atlas {
namespace util {

typedef unsigned long checksum_t;

checksum_t checksum(const int values[], size_t size);
checksum_t checksum(const long values[], size_t size);
checksum_t checksum(const float values[], size_t size);
checksum_t checksum(const double values[], size_t size);
checksum_t checksum(const checksum_t values[], size_t size);

void checksum_reset(checksum_t&);
void checksum_update(checksum_t&, const int values[], size_t size);
void checksum_update(checksum_t&, const long values[], size_t size);
void checksum_update(checksum_t&, const float values[], size_t size);
void checksum_update(checksum_t&, const double values[], size_t size);
void checksum_update(checksum_t&, const checksum_t values[], size_t size);
void checksum_update(checksum_t&, const int values[], size_t size, size_t stride);
void checksum_update(checksum_t&, const long values[], size_t size, size_t stride);
void checksum_update(checksum_t&, const float values[], size_t size, size_t stride);
void checksum_update(checksum_t&, const double values[], size_t size, size_t stride);
void checksum_update(checksum_t&, const checksum_t values[], size_t size, size_t stride);
checksum_t checksum_digest(const checksum_t&);
std::string checksum_to_hex_str(const util::checksum_t& digest);

}  // namespace util
}  // namespace atlas
