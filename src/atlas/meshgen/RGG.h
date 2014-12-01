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

#include "atlas/grids/ReducedGrid.h"
#include "atlas/Util.h"

namespace atlas {
namespace meshgen {

/// @brief Reduced Gaussian Grid
class RGG : public grids::ReducedGrid
{
public:
  RGG(): grids::ReducedGrid() {}
  RGG(const int N, const int lon[]);
  RGG(const size_t N, const long lon[]);

};

/// @brief Gaussian Grid
class GG: public grids::ReducedGrid
{
public:
  GG(int nlon, int N);
};

class RegularGrid: public grids::ReducedGrid {
public:
  RegularGrid(int nlon, int nlat);
};


class T63:   public grids::ReducedGrid { public: T63();   };
class T95:   public grids::ReducedGrid { public: T95();   };
class T159:  public grids::ReducedGrid { public: T159();  };
class T255:  public grids::ReducedGrid { public: T255();  };
class T511:  public grids::ReducedGrid { public: T511();  };
class T1279: public grids::ReducedGrid { public: T1279(); };
class T2047: public grids::ReducedGrid { public: T2047(); };
class T3999: public grids::ReducedGrid { public: T3999(); };
class T7999: public grids::ReducedGrid { public: T7999(); };


grids::ReducedGrid* new_reduced_gaussian_grid(const std::string& identifier);
grids::ReducedGrid* new_reduced_gaussian_grid(const std::vector<long>& nlon);
grids::ReducedGrid* new_regular_gaussian_grid(int nlon, int nlat);
grids::ReducedGrid* new_regular_latlon_grid(int nlon, int nlat);

typedef grids::ReducedGrid __ReducedGrid;
extern "C"
{
  __ReducedGrid* atlas__new_reduced_gaussian_grid(char* identifier);
  __ReducedGrid* atlas__new_regular_gaussian_grid(int nlon, int nlat);
  __ReducedGrid* atlas__new_regular_latlon_grid(int nlon, int nlat);
  __ReducedGrid* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat);
}

} // namespace meshgen
} // namespace atlas

#endif // RGG_h
