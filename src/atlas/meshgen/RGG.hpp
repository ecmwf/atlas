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
#include "atlas/Metadata.hpp"

namespace atlas {
  class Mesh;
namespace meshgen {

struct Region;
  
class RGG
{
public:
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


class RGGMeshGenerator
{
public:
  RGGMeshGenerator();
  Mesh* generate(const RGG& rgg);
  Mesh* operator()(const RGG& rgg){ return generate(rgg); }
private:
  void generate_region(const RGG& rgg, const std::vector<int>& parts, int mypart, Region& region);
  Mesh* generate_mesh(const RGG& rgg,const std::vector<int>& parts, const Region& region);
  std::vector<int> partition(const RGG& rgg) const;
public:
  Metadata options;
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
