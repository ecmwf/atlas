
#ifndef atlas_Checksum_hpp
#define atlas_Checksum_hpp

#include <sstream>
#include <limits>

namespace atlas {

std::string checksum(const int    values[], size_t size);
std::string checksum(const float  values[], size_t size);
std::string checksum(const double values[], size_t size);

} // namespace atlas

#endif
