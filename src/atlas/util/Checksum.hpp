
#ifndef atlas_Checksum_hpp
#define atlas_Checksum_hpp

#include <sstream>
#include <limits>

namespace atlas {

typedef unsigned long checksum_t;

checksum_t checksum(const int    values[], size_t size);
checksum_t checksum(const long   values[], size_t size);
checksum_t checksum(const float  values[], size_t size);
checksum_t checksum(const double values[], size_t size);
checksum_t checksum(const checksum_t  values[], size_t size);

std::string checksum_str(const int    values[], size_t size);
std::string checksum_str(const long   values[], size_t size);
std::string checksum_str(const float  values[], size_t size);
std::string checksum_str(const double values[], size_t size);
std::string checksum_str(const checksum_t values[], size_t size);

} // namespace atlas

#endif
