/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef RGG_hpp
#define RGG_hpp

#include <math.h>
#include <vector>

namespace atlas {
namespace meshgen {
  
/// @brief Reduced Gaussian Grid
class RGG
{
public:
  RGG() {}
  RGG(const int nlat, const int lon[]);
  RGG(const size_t nlat, const long lon[]);
  int nlat() const { return lat_.size(); }
  int nlon(int jlat) const { return lon_[jlat]; }
  int nlonmax() const { return lon_[nlat()/2]; }
  double lon(const int jlon, const int jlat) const { return 2.*M_PI/static_cast<double>(nlon(jlat))*static_cast<double>(jlon); }
  double lat(const int jlat) const { return lat_[jlat]; }
  int ngptot() const;
protected:
  std::vector<double> lat_;
  std::vector<int>    lon_;
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


} // namespace meshgen
} // namespace atlas

#endif // RGG_hpp
