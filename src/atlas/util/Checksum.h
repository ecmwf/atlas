#pragma once

#include <cstddef>

namespace atlas {
namespace util {

typedef unsigned long checksum_t;

checksum_t checksum( const int values[], size_t size );
checksum_t checksum( const long values[], size_t size );
checksum_t checksum( const float values[], size_t size );
checksum_t checksum( const double values[], size_t size );
checksum_t checksum( const checksum_t values[], size_t size );

}  // namespace util
}  // namespace atlas
