
#include <iostream>
#include "atlas/private/Checksum.h"

namespace atlas {

namespace
{
typedef unsigned long  uint64_t;
typedef unsigned int   uint32_t;
typedef unsigned short uint16_t;


// Unused private function
#if 0
uint32_t fletcher32( uint16_t const *data, size_t words )
{
  uint32_t sum1 = 0xffff, sum2 = 0xffff;

  while (words) {
          unsigned tlen = words > 359 ? 359 : words;
          words -= tlen;
          do {
                  sum2 += sum1 += *data++;
          } while (--tlen);
          sum1 = (sum1 & 0xffff) + (sum1 >> 16);
          sum2 = (sum2 & 0xffff) + (sum2 >> 16);
  }
  /* Second reduction step to reduce sums to 16 bits */
  sum1 = (sum1 & 0xffff) + (sum1 >> 16);
  sum2 = (sum2 & 0xffff) + (sum2 >> 16);
  return sum2 << 16 | sum1;
}
#endif


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
//  return fletcher32( reinterpret_cast<const uint16_t*>(data), size/sizeof(uint16_t) );
}

checksum_t checksum(const int values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(int)/sizeof(char));
}

checksum_t checksum(const long values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(long)/sizeof(char));
}

checksum_t checksum(const float values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(float)/sizeof(char));
}

checksum_t checksum(const double values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(double)/sizeof(char));
}

checksum_t checksum(const checksum_t values[], size_t size)
{
  return checksum(reinterpret_cast<const char*>(&values[0]),size*sizeof(checksum_t)/sizeof(char));
}

}// namespace atlas
