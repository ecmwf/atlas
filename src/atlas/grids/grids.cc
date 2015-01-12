/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/parser/Tokenizer.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/Builder.h"
#include "atlas/grids/grids.h"
#include "atlas/GridSpecParams.h"

using eckit::Log;
using eckit::Tokenizer;
using eckit::Translator;
using eckit::Factory;

namespace atlas {
namespace grids {

Grid* grid_from_uid(const std::string& uid)
{
  Grid* grid = NULL;
  if( Factory<Grid>::instance().exists(uid) )
  {
    grid = Grid::create( GridSpec(uid) );
  }
  else
  {
    Tokenizer tokenize(".");
    std::vector<std::string> tokens;
    tokenize(uid,tokens);
    Translator<std::string,int> to_int;
    std::string grid_type = tokens[0];
    if( grid_type == "ll" ) grid_type = LonLatGrid::grid_type_str();
    if( grid_type == "gg" ) grid_type = GaussianGrid::grid_type_str();
    if( grid_type == "rgg") grid_type = ReducedGaussianGrid::grid_type_str();

    if( grid_type == ReducedGaussianGrid::grid_type_str() )
    {
      throw eckit::BadParameter("Grid ["+uid+"] does not exist.",Here());
    }
    else if( tokens.size() > 1)
    {
      GridSpec gridspec(grid_type);

      if( tokens[1][0] == 'N' )
      {
        std::string Nstr(tokens[1],1,tokens[1].size()-1);
        int N = to_int(Nstr);
        gridspec.set("N",N);
      }
      else
      {
        std::vector<std::string> lonlat;
        Tokenizer tokenize_lonlat("x");
        tokenize_lonlat(tokens[1],lonlat);
        if( lonlat.size() > 1 )
        {
          int nlon = to_int(lonlat[0]);
          int nlat = to_int(lonlat[1]);
          gridspec.set("nlon",nlon);
          gridspec.set("nlat",nlat);
        }
      }
      grid = Grid::create(gridspec);
    }
    else
    {
      throw eckit::BadParameter("Insufficient information to construct grid "+uid+" or grid does not exist.",Here());
    }
  }
  return grid;
}


template<typename CONCRETE>
void load_grid()
{
  eckit::ConcreteBuilderT1<Grid,CONCRETE> builder("tmp");
}

void load()
{
  eckit::Log::debug() << "Loading library [atlas::grids]" << std::endl;

  // We have to touch all classes we want to register for static linking.

  load_grid<ReducedGrid>();
  load_grid<GaussianGrid>();
  load_grid<ReducedGaussianGrid>();
  load_grid<LonLatGrid>();
  load_grid<ReducedLonLatGrid>();
  load_grid<RotatedLatLon>();
  load_grid<Unstructured>();
  load_grid<PolarStereoGraphic>();

  load_grid<rgg::N16>();
  load_grid<rgg::N24>();
  load_grid<rgg::N32>();
  load_grid<rgg::N48>();
  load_grid<rgg::N64>();
  load_grid<rgg::N80>();
  load_grid<rgg::N96>();
  load_grid<rgg::N128>();
  load_grid<rgg::N160>();
  load_grid<rgg::N200>();
  load_grid<rgg::N256>();
  load_grid<rgg::N320>();
  load_grid<rgg::N400>();
  load_grid<rgg::N512>();
  load_grid<rgg::N576>();
  load_grid<rgg::N640>();
  load_grid<rgg::N800>();
  load_grid<rgg::N1024>();
  load_grid<rgg::N1280>();
  load_grid<rgg::N1600>();
  load_grid<rgg::N2000>();
  load_grid<rgg::N4000>();
  load_grid<rgg::N8000>();

  load_grid<rgg::OctahedralRGG>();

}

} // namespace grids
} // namespace atlas

extern "C"
{
  void atlas__grids__load()
  {
    atlas::grids::load();
  }
}
