/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Parameters_h
#define atlas_Parameters_h

#include <cmath>

#include "atlas/atlas_defines.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

enum { XX = 0, YY = 1, ZZ = 2 };

enum { LON = 0, LAT = 1 };

/// @brief The usual PI/180 constant
static const double DEG_TO_RAD = M_PI/180.;

//------------------------------------------------------------------------------------------------------

struct EARTH
{
  static const double RADIUS = 6371.22e+03;
};

struct Entity
{
		enum { NODES=0, FACES=1, ELEMS=2 };
};

//------------------------------------------------------------------------------------------------------

struct ElementRef
{

		ElementRef() {}

		ElementRef(int func_space_idx, int elem_idx) :f(func_space_idx),e(elem_idx)
		{
		}

		int f;
		int e;

		bool operator<(const ElementRef& other) const
		{
			if( f <  other.f ) return true;
			if( f == other.f ) return (e < other.e);
			return false;
		};
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_Parameters_h
