
#ifndef atlas_Checksum_h
#define atlas_Checksum_h

#include <sstream>
#include <limits>

namespace atlas {
namespace internals {
typedef unsigned long checksum_t;

checksum_t checksum(const int    values[], size_t size);
checksum_t checksum(const long   values[], size_t size);
checksum_t checksum(const float  values[], size_t size);
checksum_t checksum(const double values[], size_t size);
checksum_t checksum(const checksum_t  values[], size_t size);

} // namespace internals
} // namespace atlas

#endif
