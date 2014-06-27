
#include "atlas/util/Checksum.hpp"

namespace atlas {

namespace
{
typedef unsigned long uint64_t;
typedef unsigned int  uint32_t;

inline uint64_t fletcher64( const uint32_t* data, size_t count )
{
  uint64_t sum1 = 0xffffffff, sum2 = 0xffffffff;

  while (count) {
    size_t tlen = count > 92679 ? 92679 : count;
    count -= tlen;
    do {
      sum2 += sum1 += *data++;
    } while (--tlen);
    sum1 = (sum1 & 0xffffffff) + (sum1 >> 32);
    sum2 = (sum2 & 0xffffffff) + (sum2 >> 32);
  }
  /* Second reduction step to reduce sums to 32 bits */
  sum1 = (sum1 & 0xffffffff) + (sum1 >> 32);
  sum2 = (sum2 & 0xffffffff) + (sum2 >> 32);
  return sum2 << 32 | sum1;
}

}

std::string checksum(const char* data, size_t size)
{
  std::stringstream ss;
  ss << fletcher64( reinterpret_cast<const uint32_t*>(data), size/sizeof(uint32_t) );
  return ss.str();
}

std::string checksum(const int values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(int));
}

std::string checksum(const float values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(float));
}

std::string checksum(const double values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(double));
}

}// namespace atlas
