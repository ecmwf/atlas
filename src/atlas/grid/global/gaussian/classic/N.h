/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Mar 2016

#ifndef atlas_grids_global_gaussian_classic_N_h
#define atlas_grids_global_gaussian_classic_N_h

#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"

namespace atlas {
namespace grid {
namespace global {
namespace gaussian {
namespace classic {

class PointsPerLatitude : public eckit::Owned
{
public:
  typedef eckit::BuilderT0<PointsPerLatitude> builder_t;

  static std::string className()
    { return "atlas.global.gaussian.classic.PointsPerLatitude"; }

  /// @pre nlats has enough allocated memory to store the latitudes
  /// @param size of lats vector
  void assign(long nlon[], const size_t size) const;

  /// @post resizes the vector to the number of latitutes
  void assign(std::vector<long>& nlon) const;

  size_t N() const { return nlon_.size(); }

protected:

  std::vector<long> nlon_;

};

//-----------------------------------------------------------------------------

#define DECLARE_POINTS_PER_LATITUDE(NUMBER) \
  class N##NUMBER : public PointsPerLatitude { public: N##NUMBER(); };

#define LIST(...) __VA_ARGS__
#define DEFINE_POINTS_PER_LATITUDE(NUMBER,NLON) \
  eckit::ConcreteBuilderT0<PointsPerLatitude,N##NUMBER> builder_N##NUMBER(#NUMBER); \
  \
  N##NUMBER::N##NUMBER()\
  {\
    size_t N = NUMBER;\
    long nlon[] = {NLON} ;\
    nlon_.assign(nlon,nlon+N);\
  }\

//-----------------------------------------------------------------------------

DECLARE_POINTS_PER_LATITUDE(16);
DECLARE_POINTS_PER_LATITUDE(24);
DECLARE_POINTS_PER_LATITUDE(32);
DECLARE_POINTS_PER_LATITUDE(48);
DECLARE_POINTS_PER_LATITUDE(64);
DECLARE_POINTS_PER_LATITUDE(80);
DECLARE_POINTS_PER_LATITUDE(96);
DECLARE_POINTS_PER_LATITUDE(128);
DECLARE_POINTS_PER_LATITUDE(160);
DECLARE_POINTS_PER_LATITUDE(200);
DECLARE_POINTS_PER_LATITUDE(256);
DECLARE_POINTS_PER_LATITUDE(320);
DECLARE_POINTS_PER_LATITUDE(400);
DECLARE_POINTS_PER_LATITUDE(512);
DECLARE_POINTS_PER_LATITUDE(576);
DECLARE_POINTS_PER_LATITUDE(640);
DECLARE_POINTS_PER_LATITUDE(800);
DECLARE_POINTS_PER_LATITUDE(1024);
DECLARE_POINTS_PER_LATITUDE(1280);
DECLARE_POINTS_PER_LATITUDE(1600);
DECLARE_POINTS_PER_LATITUDE(2000);
DECLARE_POINTS_PER_LATITUDE(4000);
DECLARE_POINTS_PER_LATITUDE(8000);

#undef DECLARE_POINTS_PER_LATITUDE

//-----------------------------------------------------------------------------

} // namespace classic
} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas

#include "atlas/grid/global/gaussian/Gaussian.h"
#include "atlas/grid/global/gaussian/latitudes/Latitudes.h"

namespace atlas {
namespace grid {
namespace global {
namespace rgg {

//------------------------------------------------------------------------------------------------------

#define DEFINE_GRID(CLASS)\
class CLASS : public ReducedGaussianGrid { \
public:\
\
  static std::string grid_type_str() { return "rgg."+std::string(#CLASS); } \
  CLASS() \
  {\
    construct();\
    ReducedGrid::N_ = nlat()/2;\
    set_typeinfo();\
  }\
  CLASS(Grid::ARG1 arg1)\
  {\
    construct();\
    ReducedGrid::N_ = nlat()/2;\
    set_typeinfo();\
  }\
  void construct();\
  void set_typeinfo() { \
    shortName_ = std::string(#CLASS); \
    grid_type_ = ReducedGaussianGrid::grid_type_str(); \
  }\
  static std::string className() { return "atlas.grids.rgg."+std::string(#CLASS); }\
};

// list of reduced_gg grid definitions

DEFINE_GRID(N16);
DEFINE_GRID(N24);
DEFINE_GRID(N32);
DEFINE_GRID(N48);
DEFINE_GRID(N64);
DEFINE_GRID(N80);
DEFINE_GRID(N96);
DEFINE_GRID(N128);
DEFINE_GRID(N160);
DEFINE_GRID(N200);
DEFINE_GRID(N256);
DEFINE_GRID(N320);
DEFINE_GRID(N400);
DEFINE_GRID(N512);
DEFINE_GRID(N576);
DEFINE_GRID(N640);
DEFINE_GRID(N800);
DEFINE_GRID(N1024);
DEFINE_GRID(N1280);
DEFINE_GRID(N1600);
DEFINE_GRID(N2000);
DEFINE_GRID(N4000);
DEFINE_GRID(N8000);

#undef DEFINE_GRID

//------------------------------------------------------------------------------------------------------

} // namespace rgg
} // namespace global
} // namespace grid
} // namespace atlas

#endif // atlas_grids_global_gaussian_classic_N_h
