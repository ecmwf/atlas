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
/// @date Nov 2014

#ifndef atlas_grids_gausslat_gausslat_h
#define atlas_grids_gausslat_gausslat_h

#include "eckit/memory/Owned.h"
#include "eckit/memory/Builder.h"

namespace atlas {
namespace grids {
namespace gausslat {

//------------------------------------------------------------------------------------------------------

class GaussianLatitudes : public eckit::Owned
{
public:
  typedef eckit::BuilderT0<GaussianLatitudes> builder_t;

  static std::string className() { return "GaussianLatitudes"; }

  /// @pre nlats has enough allocated memory to store the latitudes
  /// @param size of lats vector
  void assign(double lats[], const size_t size) const;

  /// @post resizes the vector to the number of latitutes
  void assign(std::vector<double>& lats) const;

  size_t N() const { return lats_.size(); }

protected:

  std::vector<double> lats_;

};

#define DECLARE_GAUSSIAN_LATITUDES(NUMBER) \
  class N##NUMBER : public GaussianLatitudes { public: N##NUMBER(); };

#define LIST(...) __VA_ARGS__
#define DEFINE_GAUSSIAN_LATITUDES(NUMBER,LATS) \
  eckit::ConcreteBuilderT0<GaussianLatitudes,N##NUMBER> builder_N##NUMBER(#NUMBER); \
  \
  N##NUMBER::N##NUMBER()\
  {\
    size_t N = NUMBER;\
    double lat[] = {LATS} ;\
    lats_.assign(lat,lat+N);\
  }\

DECLARE_GAUSSIAN_LATITUDES(16);
DECLARE_GAUSSIAN_LATITUDES(24);
DECLARE_GAUSSIAN_LATITUDES(32);
DECLARE_GAUSSIAN_LATITUDES(48);
DECLARE_GAUSSIAN_LATITUDES(64);
DECLARE_GAUSSIAN_LATITUDES(80);
DECLARE_GAUSSIAN_LATITUDES(96);
DECLARE_GAUSSIAN_LATITUDES(128);
DECLARE_GAUSSIAN_LATITUDES(160);
DECLARE_GAUSSIAN_LATITUDES(200);
DECLARE_GAUSSIAN_LATITUDES(256);
DECLARE_GAUSSIAN_LATITUDES(320);
DECLARE_GAUSSIAN_LATITUDES(400);
DECLARE_GAUSSIAN_LATITUDES(512);
DECLARE_GAUSSIAN_LATITUDES(576);
DECLARE_GAUSSIAN_LATITUDES(640);
DECLARE_GAUSSIAN_LATITUDES(800);
DECLARE_GAUSSIAN_LATITUDES(1024);
DECLARE_GAUSSIAN_LATITUDES(1280);
DECLARE_GAUSSIAN_LATITUDES(1600);
DECLARE_GAUSSIAN_LATITUDES(2000);
DECLARE_GAUSSIAN_LATITUDES(4000);
DECLARE_GAUSSIAN_LATITUDES(8000);

#undef DECLARE_GAUSSIAN_LATITUDES

//------------------------------------------------------------------------------------------------------

} // namespace gausslat
} // namespace grids
} // namespace atlas

#endif // atlas_grids_gausslat_gausslat_h
