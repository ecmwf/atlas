/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Nov 2014

#ifndef atlas_grids_reduced_gg_reduced_gg_h
#define atlas_grids_reduced_gg_reduced_gg_h

#include "eckit/memory/Builder.h"
#include "atlas/grids/ReducedGaussianGrid.h"
#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/grids/rgg/OctahedralRGG.h"

namespace atlas {
namespace grids {
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
} // namespace grids
} // namespace atlas

#endif // atlas_grids_reduced_gg_reduced_gg_h
