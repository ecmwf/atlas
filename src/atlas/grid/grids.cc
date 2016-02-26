/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <regex.h>
#include "eckit/parser/Tokenizer.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/Builder.h"
#include "atlas/grid/grids.h"
#include "atlas/util/runtime/Log.h"

using eckit::Tokenizer;
using eckit::Translator;
using eckit::Factory;

namespace atlas {
namespace grid {
  
size_t regex_count_parens(const std::string& string)
{
    size_t out = 0;
    bool last_was_backslash = 0;
    for(const char *step=string.c_str(); *step !='\0'; step++){
        if (*step == '\\' && !last_was_backslash){
            last_was_backslash = true;
            continue;
        }
        if (*step == ')' && !last_was_backslash)
            out++;
        last_was_backslash = false;
    }
    return out;
}

int regex_match_impl( const std::string& string, 
                      const std::string& regex, 
                      std::vector<std::string>& substr, 
                      bool use_substr, 
                      bool use_case)
{
    regex_t re;
    size_t matchcount = 0;
    if (use_substr)
      matchcount = regex_count_parens(regex);
    regmatch_t result[matchcount+1];
    int compiled_ok = !regcomp(&re, regex.c_str(), REG_EXTENDED
                                            + (use_case   ? 0 : REG_ICASE)
                                            + (use_substr ? 0 : REG_NOSUB) );

    if( !compiled_ok ) 
      printf("This regular expression didn't compile: \"%s\"", regex.c_str());

    int found = !regexec(&re, string.c_str(), matchcount+1, result, 0);
    if (found && use_substr){
        substr.resize(matchcount);
        //match zero is the whole string; ignore it.
        for (size_t i=0; i< matchcount; i++){
            if (result[i+1].rm_eo > 0){ //GNU peculiarity: match-to-empty marked with -1.
                size_t length_of_match = result[i+1].rm_eo - result[i+1].rm_so;
                substr[i] = std::string(&string[result[i+1].rm_so],length_of_match);
            }
        }
    }
    regfree(&re);
    return found;
}

class Regex
{
public:
  Regex(const std::string& regex, bool use_case=true) : regex_(regex), use_case_(use_case) {}
  bool match(const std::string& string)
  {
    std::vector<std::string> substr;
    return regex_match_impl(string,regex_,substr,false,use_case_);
  }
  bool match(const std::string& string, std::vector<std::string>& substr)
  {
    return regex_match_impl(string,regex_,substr,true,use_case_);
  }
private:
  std::string regex_;
  bool use_case_;
};

//=====================================================

Grid* grid_from_uid(const std::string& uid)
{
  if( Factory<Grid>::instance().exists(uid) )
  {
    return Grid::create(util::Config("grid_type", uid));
  }
  else
  {
    Regex classical_reduced_gaussian_grid  ("^N([0-9]+)$");
    Regex octahedral_reduced_gaussian_grid ("^O([0-9]+)$");
    Regex regular_gaussian_grid            ("^F([0-9]+)$");
    
    util::Config gridparams;
    Translator<std::string,int> to_int;
    std::vector<std::string> matches;
    if( classical_reduced_gaussian_grid.match(uid) )
    {
      throw eckit::BadParameter("Grid ["+uid+"] does not exist.",Here());
    }
    else if( octahedral_reduced_gaussian_grid.match(uid,matches) )
    {
      int N = to_int(matches[0]);
      gridparams.set("grid_type", OctahedralRGG::grid_type_str());
      gridparams.set("N",N);
      return Grid::create( gridparams );
    }
    else if( regular_gaussian_grid.match(uid,matches) )
    {
      int N = to_int(matches[0]);
      gridparams.set("grid_type", GaussianGrid::grid_type_str());
      gridparams.set("N",N);
      return Grid::create( gridparams );
    }
    else
    {
      Tokenizer tokenize(".");
      std::vector<std::string> tokens;
      tokenize(uid,tokens);
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
        gridparams.set("grid_type",grid_type);
        if( tokens[1][0] == 'N' )
        {
          std::string Nstr(tokens[1],1,tokens[1].size()-1);
          int N = to_int(Nstr);
          gridparams.set("N",N);
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
            gridparams.set("nlon",nlon);
            gridparams.set("nlat",nlat);
            if( nlat%2 == 1 ) gridparams.set("poles",true);
          }
        }
      }
      return Grid::create( gridparams );
    }
  }
  throw eckit::BadParameter("Insufficient information to construct grid "+uid+" or grid does not exist.",Here());
  return 0;
}


template<typename CONCRETE>
void load_grid()
{
  eckit::ConcreteBuilderT1<Grid,CONCRETE> builder("tmp");
}

void load()
{
  Log::debug(2) << "Loading library [atlas::grid]" << std::endl;

  // We have to touch all classes we want to register for static linking.

  load_grid<ReducedGrid>();
  load_grid<GaussianGrid>();
  load_grid<ReducedGaussianGrid>();
  load_grid<LonLatGrid>();
  load_grid<ReducedLonLatGrid>();
  load_grid<Unstructured>();

  load_grid<predefined::rgg::N16>();
  load_grid<predefined::rgg::N24>();
  load_grid<predefined::rgg::N32>();
  load_grid<predefined::rgg::N48>();
  load_grid<predefined::rgg::N64>();
  load_grid<predefined::rgg::N80>();
  load_grid<predefined::rgg::N96>();
  load_grid<predefined::rgg::N128>();
  load_grid<predefined::rgg::N160>();
  load_grid<predefined::rgg::N200>();
  load_grid<predefined::rgg::N256>();
  load_grid<predefined::rgg::N320>();
  load_grid<predefined::rgg::N400>();
  load_grid<predefined::rgg::N512>();
  load_grid<predefined::rgg::N576>();
  load_grid<predefined::rgg::N640>();
  load_grid<predefined::rgg::N800>();
  load_grid<predefined::rgg::N1024>();
  load_grid<predefined::rgg::N1280>();
  load_grid<predefined::rgg::N1600>();
  load_grid<predefined::rgg::N2000>();
  load_grid<predefined::rgg::N4000>();
  load_grid<predefined::rgg::N8000>();

  load_grid<OctahedralRGG>();

}


//extern "C"
//{

  ReducedGrid* new_reduced_grid( const std::string& identifier )
  {
    return ReducedGrid::create(identifier);
  }

  ReducedGrid* new_reduced_gaussian_grid( const std::vector<long>& nlon )
  {
    return new ReducedGaussianGrid(nlon.size(),nlon.data());
  }

  ReducedGrid* new_lonlat_grid( int nlon, int nlat )
  {
    return new LonLatGrid(static_cast<size_t>(nlon),static_cast<size_t>(nlat));
  }

  ReducedGrid* new_gaussian_grid( int N )
  {
    return new GaussianGrid(N);
  }


  ReducedGrid* atlas__new_reduced_grid(char* identifier)
  {
    return new_reduced_grid( std::string(identifier) );
  }

  ReducedGrid* atlas__new_gaussian_grid ( int N )
  {
    return new_gaussian_grid( N );
  }

  ReducedGrid* atlas__new_lonlat_grid(int nlon, int nlat)
  {
    return new_lonlat_grid( nlon, nlat );
  }

  ReducedGrid* atlas__new_reduced_gaussian_grid(int nlon[], int nlat)
  {
    std::vector<long> nlon_vector;
    nlon_vector.assign(nlon,nlon+nlat);
    return new_reduced_gaussian_grid(nlon_vector);
  }

  void atlas__grids__load()
  {
    atlas::grid::load();
  }

  void atlas__ReducedGrid__delete(ReducedGrid* This)
  {
    delete This;
  }

//}

} // namespace grid
} // namespace atlas

