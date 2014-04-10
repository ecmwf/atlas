// (C) Copyright 1996-2014 ECMWF.

#ifndef RGG_hpp
#define RGG_hpp

#include <math.h>
#include <vector>
#include "atlas/Metadata.hpp"

namespace atlas {
  class Mesh;
namespace meshgen {

// For Hemisphere!, jlat=0 at equator
class RGG
{
public:
  int nlat() const { return lat_.size(); }
  int nlon(int jlat) const { return lon_[jlat]; }
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

} // namespace meshgen
} // namespace atlas

#endif // RGG_hpp
