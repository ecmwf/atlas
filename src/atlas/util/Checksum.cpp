
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

checksum_t checksum(const char* data, size_t size)
{
  return fletcher64( reinterpret_cast<const uint32_t*>(data), size/sizeof(uint32_t) );
}

checksum_t checksum(const int values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(int));
}

checksum_t checksum(const long values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(long));
}

checksum_t checksum(const float values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(float));
}

checksum_t checksum(const double values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(double));
}

checksum_t checksum(const checksum_t values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(checksum_t));
}


std::string checksum_str(const int values[], size_t size)
{
  std::stringstream ss; ss << checksum(values,size);
  return ss.str();
}

std::string checksum_str(const long values[], size_t size)
{
  std::stringstream ss; ss << checksum(values,size);
  return ss.str();
}

std::string checksum_str(const float values[], size_t size)
{
  std::stringstream ss; ss << checksum(values,size);
  return ss.str();
}


std::string checksum_str(const double values[], size_t size)
{
  std::stringstream ss; ss << checksum(values,size);
  return ss.str();
}

std::string checksum_str(const checksum_t values[], size_t size)
{
  std::stringstream ss; ss << checksum(values,size);
  return ss.str();
}


}// namespace atlas
