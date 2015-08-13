/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an size_tergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_DataType_h
#define atlas_DataType_h

#include <string>
#include <sstream>

#include "eckit/exception/Exceptions.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

struct DataType {
  static const long KIND_INT32  = -4;
  static const long KIND_INT64  = -8;
  static const long KIND_REAL32 =  4;
  static const long KIND_REAL64 =  8;
  static std::string int32();
  static std::string int64();
  static std::string real32();
  static std::string real64();

  template< typename DATATYPE > static long kind();
  template< typename DATATYPE > static long kind(const DATATYPE&);
  
  template< typename DATATYPE > static std::string datatype();
  template< typename DATATYPE > static std::string datatype(const DATATYPE);

  static long datatype_to_kind(const std::string&);
  static std::string kind_to_datatype(const long&);
  static bool kind_valid(const long&);
};
inline std::string DataType::int32()  { return "int32";  }
inline std::string DataType::int64()  { return "int64";  }
inline std::string DataType::real32() { return "real32"; }
inline std::string DataType::real64() { return "real64"; }
template<> inline std::string DataType::datatype<int>()    { return int32();  }
template<> inline std::string DataType::datatype<long>()   { return int64();  }
template<> inline std::string DataType::datatype<float>()  { return real32(); }
template<> inline std::string DataType::datatype<double>() { return real64(); }
template<> inline std::string DataType::datatype(const int&)    { return int32();  }
template<> inline std::string DataType::datatype(const long&)   { return int64();  }
template<> inline std::string DataType::datatype(const float&)  { return real32(); }
template<> inline std::string DataType::datatype(const double&) { return real64(); }
template<> inline long DataType::kind<int>()    { return KIND_INT32;    }
template<> inline long DataType::kind<long>()   { return KIND_INT64;    }
template<> inline long DataType::kind<float>()  { return KIND_REAL32;   }
template<> inline long DataType::kind<double>() { return KIND_REAL64;   }
template<> inline long DataType::kind(const int&)    { return KIND_INT32;   }
template<> inline long DataType::kind(const long&)   { return KIND_INT64;   }
template<> inline long DataType::kind(const float&)  { return KIND_REAL32;   }
template<> inline long DataType::kind(const double&) { return KIND_REAL64;   }

inline long DataType::datatype_to_kind(const std::string& datatype)
{
  if      ( datatype == "int32"  ) return KIND_INT32;
  else if ( datatype == "int64"  ) return KIND_INT64;
  else if ( datatype == "real32" ) return KIND_REAL32;
  else if ( datatype == "real64" ) return KIND_REAL64;
  else throw eckit::Exception("datatype "+datatype+" not recognised.");
  return 0;
}
inline std::string DataType::kind_to_datatype(const long &kind)
{
  switch( kind )
  {
    case KIND_INT32: return int32();
    case KIND_INT64: return int64();
    case KIND_REAL32: return real32();
    case KIND_REAL64: return real64();
    default:
    std::stringstream msg;
    msg << "kind "<<kind<<" not recognised.";
    throw eckit::Exception(msg.str());
  }
}
inline bool DataType::kind_valid(const long &kind)
{
  switch( kind )
  {
    case KIND_INT32:
    case KIND_INT64:
    case KIND_REAL32:
    case KIND_REAL64:
      return true;
    default:
      return false;
  }
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
