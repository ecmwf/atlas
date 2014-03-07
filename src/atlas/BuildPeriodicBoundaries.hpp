// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#ifndef BuildPeriodicBoundaries_hpp
#define BuildPeriodicBoundaries_hpp
#include <string>
namespace atlas {
class Mesh;

void build_periodic_boundaries( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void atlas__build_periodic_boundaries (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // BuildPeriodicBoundaries_hpp
