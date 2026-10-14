/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Checksum.h"

namespace atlas {
namespace util {
namespace {  // anonymous

template <typename Value>
const char* as_bytes(const Value values[], size_t size) {
    return size ? reinterpret_cast<const char*>(values) : nullptr;
}

template <typename Value>
size_t byte_size(size_t size) {
    return size * sizeof(Value) / sizeof(char);
}

union Fletcher16 {
  checksum_t checksum;
  uint16_t c[2];
};

static void fletcher16_reset(Fletcher16& f) {
  f.c[0] = 0;
  f.c[1] = 0;
}

static void fletcher16_update(Fletcher16& f, const uint8_t* data, size_t size) {
  uint32_t c0 = f.c[0];
  uint32_t c1 = f.c[1];
  while(size > 0) {
    size_t blocklen = size;
    if (blocklen > 5802) {
      blocklen = 5802;
    }
    size -= blocklen;
    do {
      c0 = c0 + *data++;
      c1 = c1 + c0;
    } while (--blocklen);
    c0 = c0 % 0xff;
    c1 = c1 % 0xff;
  }
  f.c[0] = c0;
  f.c[1] = c1;
}

static void fletcher16_update(Fletcher16& f, const uint8_t* data, size_t size, size_t stride, size_t chunk) {
  if (stride == 1 && chunk == 1) {
    fletcher16_update(f, data, size);
    return;
  }

  uint32_t c0 = f.c[0];
  uint32_t c1 = f.c[1];
  const size_t max_blocklen = std::max(size_t{1}, size_t{5802} / chunk);
  while (size > 0) {
    size_t blocklen = size;
    if (blocklen > max_blocklen) {
      blocklen = max_blocklen;
    }
    size -= blocklen;
    do {
      for (size_t byte = 0; byte < chunk; ++byte) {
        c0 = c0 + data[byte];
        c1 = c1 + c0;
      }
      data += stride;
    } while (--blocklen);
    c0 = c0 % 0xff;
    c1 = c1 % 0xff;
  }
  f.c[0] = c0;
  f.c[1] = c1;
}

static uint16_t fletcher16_finish(const Fletcher16& f) {
  uint32_t c0 = f.c[0];
  uint32_t c1 = f.c[1];
  return (c1 << 8 | c0);
}

static uint16_t fletcher16(const uint8_t* data, size_t size) {
  Fletcher16 checksum;
  fletcher16_reset(checksum);
  fletcher16_update(checksum, data, size);
  return fletcher16_finish(checksum);
}

static void checksum_update(checksum_t& checksum, const char* data, size_t bytes) {
  fletcher16_update(reinterpret_cast<Fletcher16&>(checksum), reinterpret_cast<const uint8_t*>(data), bytes);
}

template <typename Value>
static void checksum_update_strided(checksum_t& checksum, const Value values[], size_t size, size_t stride) {
  auto& f = reinterpret_cast<Fletcher16&>(checksum);
  const auto* data = reinterpret_cast<const uint8_t*>(values);
  const size_t value_bytes = sizeof(Value);
  const size_t stride_bytes = stride * value_bytes;
  fletcher16_update(f, data, size, stride_bytes, value_bytes);
}

static checksum_t checksum(const char* data, size_t size) {
    return fletcher16(reinterpret_cast<const uint8_t*>(data), size / sizeof(uint8_t));
}

}  // namespace

void checksum_reset(checksum_t& checksum) {
  fletcher16_reset(reinterpret_cast<Fletcher16&>(checksum));
}

checksum_t checksum(const int values[], size_t size) {
  return checksum(as_bytes(values, size), byte_size<int>(size));
}

checksum_t checksum(const long values[], size_t size) {
  return checksum(as_bytes(values, size), byte_size<long>(size));
}

checksum_t checksum(const float values[], size_t size) {
  return checksum(as_bytes(values, size), byte_size<float>(size));
}

checksum_t checksum(const double values[], size_t size) {
  return checksum(as_bytes(values, size), byte_size<double>(size));
}

checksum_t checksum(const checksum_t values[], size_t size) {
  return checksum(as_bytes(values, size), byte_size<checksum_t>(size));
}

void checksum_update(checksum_t& c, const int values[], size_t size) {
  checksum_update(c, as_bytes(values, size), byte_size<int>(size));
}

void checksum_update(checksum_t& c, const long values[], size_t size) {
  checksum_update(c, as_bytes(values, size), byte_size<long>(size));
}

void checksum_update(checksum_t& c, const float values[], size_t size) {
  checksum_update(c, as_bytes(values, size), byte_size<float>(size));
}

void checksum_update(checksum_t& c, const double values[], size_t size) {
  checksum_update(c, as_bytes(values, size), byte_size<double>(size));
}

void checksum_update(checksum_t& c, const checksum_t values[], size_t size) {
  checksum_update(c, as_bytes(values, size), byte_size<checksum_t>(size));
}

void checksum_update(checksum_t& c, const int values[], size_t size, size_t stride) {
  checksum_update_strided(c, values, size, stride);
}

void checksum_update(checksum_t& c, const long values[], size_t size, size_t stride) {
  checksum_update_strided(c, values, size, stride);
}

void checksum_update(checksum_t& c, const float values[], size_t size, size_t stride) {
  checksum_update_strided(c, values, size, stride);
}

void checksum_update(checksum_t& c, const double values[], size_t size, size_t stride) {
  checksum_update_strided(c, values, size, stride);
}

void checksum_update(checksum_t& c, const checksum_t values[], size_t size, size_t stride) {
  checksum_update_strided(c, values, size, stride);
}

checksum_t checksum_digest(const checksum_t& checksum) {
    return fletcher16_finish(reinterpret_cast<const Fletcher16&>(checksum));
}

namespace {
inline std::string to_hex_str(util::checksum_t digest) {
  // Configure the stream for hex, lowercase, and zero padding

  std::stringstream ss;
  ss << std::hex
      << std::nouppercase
      << std::setw(4)
      << std::setfill('0')
      << digest;
  return ss.str();
}
}

std::string checksum_to_hex_str(const util::checksum_t& digest) {
  return to_hex_str(digest);
}

}  // namespace util
}  // namespace atlas
