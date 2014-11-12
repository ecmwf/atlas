/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef RGG_h
#define RGG_h

#include "atlas/ReducedGrid.h"
#include "atlas/Util.h"

namespace atlas {
namespace meshgen {

/// @brief Reduced Gaussian Grid
class RGG : public ReducedGrid
{
public:
  RGG(): ReducedGrid() {}
  RGG(const int N, const int lon[]);
  RGG(const size_t N, const long lon[]);

};

/// @brief Gaussian Grid
class GG: public ReducedGrid
{
public:
  GG(int nlon, int N);
};

class RegularGrid: public ReducedGrid {
public:
  RegularGrid(int nlon, int nlat);
};


class T63:   public ReducedGrid { public: T63();   };
class T95:   public ReducedGrid { public: T95();   };
class T159:  public ReducedGrid { public: T159();  };
class T255:  public ReducedGrid { public: T255();  };
class T511:  public ReducedGrid { public: T511();  };
class T1279: public ReducedGrid { public: T1279(); };
class T2047: public ReducedGrid { public: T2047(); };
class T3999: public ReducedGrid { public: T3999(); };
class T7999: public ReducedGrid { public: T7999(); };


ReducedGrid* new_reduced_gaussian_grid(const std::string& identifier);
ReducedGrid* new_reduced_gaussian_grid(const std::vector<long>& nlon);
ReducedGrid* new_regular_gaussian_grid(int nlon, int nlat);
ReducedGrid* new_regular_latlon_grid(int nlon, int nlat);

extern "C"
{
  ReducedGrid* atlas__new_reduced_gaussian_grid(char* identifier);
  ReducedGrid* atlas__new_regular_gaussian_grid(int nlon, int nlat);
  ReducedGrid* atlas__new_regular_latlon_grid(int nlon, int nlat);
  ReducedGrid* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat);
}

} // namespace meshgen
} // namespace atlas

#endif // RGG_h
