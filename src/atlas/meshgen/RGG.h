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

namespace atlas {
namespace meshgen {

enum AngleUnit{ DEG=0, RAD=1 };

void colat_to_lat_hemisphere(const int N, const double colat[], double lats[], const AngleUnit unit);

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[]);

void predict_gaussian_latitudes_hemisphere(const int N, double lat[]);


/// @brief Reduced Gaussian Grid
class RGG : public ReducedGrid
{
public:
  RGG(): ReducedGrid() {}
  RGG(const int nlat, const int lon[]);
  RGG(const size_t nlat, const long lon[]);

protected:
  void setup_colat_hemisphere(const int N, const int lon[], const double colat[], const AngleUnit);
  void setup_lat_hemisphere(const int N, const int lon[], const double lat[], const AngleUnit);
};

/// @brief Gaussian Grid
class GG: public RGG
{
public:
  GG(int nlon, int nlat);
};

class RegularGrid: public meshgen::RGG {
public:
  RegularGrid(int nlon, int nlat);
};


class T63:   public RGG { public: T63();   };
class T95:   public RGG { public: T95();   };
class T159:  public RGG { public: T159();  };
class T255:  public RGG { public: T255();  };
class T511:  public RGG { public: T511();  };
class T1279: public RGG { public: T1279(); };
class T2047: public RGG { public: T2047(); };
class T3999: public RGG { public: T3999(); };
class T7999: public RGG { public: T7999(); };


RGG* new_reduced_gaussian_grid(const std::string& identifier);
RGG* new_reduced_gaussian_grid(const std::vector<long>& nlon);
RGG* new_regular_gaussian_grid(int nlon, int nlat);
RGG* new_regular_latlon_grid(int nlon, int nlat);

extern "C"
{
  RGG* atlas__new_reduced_gaussian_grid(char* identifier);
  RGG* atlas__new_regular_gaussian_grid(int nlon, int nlat);
  RGG* atlas__new_regular_latlon_grid(int nlon, int nlat);
  RGG* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat);
  int  atlas__RGG__nlat(RGG* This);
  void atlas__RGG__nlon(RGG* This, const int* &nlon, int &size);
  int atlas__RGG__ngptot(RGG* This);
  double atlas__RGG__lon(RGG* This,int jlon,int jlat);
  double atlas__RGG__lat(RGG* This,int jlat);
  void atlas__RGG__lats(RGG* This, const double* &lats, int &size);
}

} // namespace meshgen
} // namespace atlas

#endif // RGG_h
