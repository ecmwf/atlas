/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_Functions_h
#define atlas_util_Functions_h

namespace atlas {
namespace util {

inline int microdeg( const double& deg )
{
  return static_cast<int>(deg*1.e6 + 0.5);
}

} // namespace util
} // namespace atlas

#endif
