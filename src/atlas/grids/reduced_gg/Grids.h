/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Nov 2014

#ifndef atlas_grids_reduced_gg_Grids_h
#define atlas_grids_reduced_gg_Grids_h

#include <eckit/memory/Builder.h>
#include <eckit/value/Params.h>
#include "atlas/ReducedGrid.h"

namespace atlas {
namespace grids {
namespace reduced_gg {


/// To get the grid:
///
///   SharedPtr< Grid > grid( Factory< Grid >::instance().get( "N80" ).create( eckit::Params() ) );
/// or
///   SharedPtr< ReducedGrid > grid( Factory< ReducedGrid >::instance().get( "N80" ).create() );


///@todo these classes will inherit from ReducedGaussianGrid when it's ready
#define DEFINE_GRID(CLASS)\
class CLASS : public ReducedGrid { \
public:\
\
  CLASS() { construct(); set_typeinfo(); }\
  CLASS(Grid::ARG1 arg1) { construct();  set_typeinfo(); }\
  void construct();\
  void set_typeinfo() { \
    uid_ = "reduced_gg."+std::string(#CLASS); \
    hash_ = uid_; \
    grid_type_ = "reduced_gg"; \
  }\
  static std::string className() { return "atlas.grids.reduced_gg."+std::string(#CLASS); }\
};\
register_BuilderT1(Grid,       CLASS,"reduced_gg."+std::string(#CLASS));\
register_BuilderT0(ReducedGrid,CLASS,"reduced_gg."+std::string(#CLASS));

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
} // namespace reduced_gg
} // namespace grids
} // namespace atlas
#endif // atlas_grids_reduced_gg_Grids_h
