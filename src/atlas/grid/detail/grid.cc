/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid.h"

#include <regex.h>
#include <iomanip>

#include "eckit/parser/Tokenizer.h"
#include "eckit/utils/Translator.h"

#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/grid/detail/grid/CustomStructured.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/grid/regular/RegularLonLat.h"
#include "atlas/grid/detail/grid/regular/ShiftedLon.h"
#include "atlas/grid/detail/grid/regular/ShiftedLat.h"
#include "atlas/grid/detail/grid/regular/ShiftedLonLat.h"
#include "atlas/grid/detail/grid/regular/RegularGaussian.h"
#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"
#include "atlas/grid/detail/grid/reduced/ReducedLonLat.h"
#include "atlas/grid/detail/grid/reduced/ClassicGaussian.h"
#include "atlas/grid/detail/grid/reduced/OctahedralGaussian.h"
#include "atlas/grid/detail/grid/Unstructured.h"
#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {

namespace {

size_t regex_count_parens(const std::string& string) {
    size_t out = 0;
    bool last_was_backslash = 0;
    for(const char *step=string.c_str(); *step !='\0'; step++) {
        if (*step == '\\' && !last_was_backslash) {
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
                      bool use_case) {
    regex_t re;
    size_t matchcount = 0;
    if (use_substr)
        matchcount = regex_count_parens(regex);
    regmatch_t result[matchcount+1];
    int compiled_ok = !regcomp(&re, regex.c_str(), REG_EXTENDED
                               + (use_case   ? 0 : REG_ICASE)
                               + (use_substr ? 0 : REG_NOSUB) );

    if( !compiled_ok )
        Log::error() << "This regular expression didn't compile: \"" << regex << "\"" << std::endl;

    ASSERT(compiled_ok);

    int found = !regexec(&re, string.c_str(), matchcount+1, result, 0);
    if (found && use_substr) {
        substr.resize(matchcount);
        //match zero is the whole string; ignore it.
        for (size_t i=0; i< matchcount; i++) {
            if (result[i+1].rm_eo > 0) {
                //GNU peculiarity: match-to-empty marked with -1.
                size_t length_of_match = result[i+1].rm_eo - result[i+1].rm_so;
                substr[i] = std::string(&string[result[i+1].rm_so],length_of_match);
            }
        }
    }
    regfree(&re);
    return found;
}

class Regex {
  public:
    Regex(const std::string& regex, bool use_case=true) :
        regex_(regex),
        use_case_(use_case) {
    }
    bool match(const std::string& string) const {
        std::vector<std::string> substr;
        return regex_match_impl(string,regex_,substr,false,use_case_);
    }
    bool match(const std::string& string, std::vector<std::string>& substr) const {
        return regex_match_impl(string,regex_,substr,true,use_case_);
    }
    operator std::string() const { return regex_; }
  private:
    std::string regex_;
    bool use_case_;
};

static eckit::Mutex *local_mutex = 0;
static GridCreator::Registry *m = 0;


static pthread_once_t once = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m = new GridCreator::Registry();
}

static eckit::Translator<std::string,int> to_int;

}  // anonymous namespace

//---------------------------------------------------------------------------------------------------------------------

  const GridCreator::Registry& GridCreator::registry() {
    return *m;
  }

  GridCreator::GridCreator( const std::string& name ) :
    name_(name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    ASSERT(registry().find(name_) == registry().end());
    (*m)[name_] = this;
  }

  GridCreator::~GridCreator() {
      pthread_once(&once, init);

      eckit::AutoLock<eckit::Mutex> lock(local_mutex);
      ASSERT(registry().find(name_) != registry().end());
      (*m).erase(name_);
  }

  const Grid::grid_t* GridCreator::create( const Grid::Config& config ) const {
    return Grid::grid_t::create( config );
  }


  bool GridCreator::match(const std::string& string, std::vector<std::string>& matches) const {
    return Regex(name_).match(string,matches);
  }

  std::string GridCreator::name() const { return name_; }

  std::ostream& operator<<( std::ostream& os, const GridCreator& g ) {
    g.print(os);
    return os;
  }

namespace {

//---------------------------------------------------------------------------------------------------------------------

static class classic_gaussian : public GridCreator {

public:

  classic_gaussian() : GridCreator("^[Nn]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "N<gauss>" << "Classic Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      Grid::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::reduced::ClassicGaussian::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} classic_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class octahedral_gaussian : public GridCreator {

public:

  octahedral_gaussian() : GridCreator("^[Oo]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::reduced::OctahedralGaussian::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} octahedral_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_gaussian : public GridCreator {

public:

  regular_gaussian() : GridCreator("^[Ff]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "F<gauss>" << "Regular Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::RegularGaussian::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} regular_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_lonlat : public GridCreator {

public:

  regular_lonlat() : GridCreator("^[Ll]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "L<gauss>" << "Regular longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::RegularLonLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} regular_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lonlat : public GridCreator {

public:

  shifted_lonlat() : GridCreator("^[Ss]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "S<gauss>" << "Shifted longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLonLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lon : public GridCreator {

public:

  shifted_lon() : GridCreator("^[Ss][Ll][Oo][Nn]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slon<gauss>" << "Shifted longitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLon::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lon_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lat : public GridCreator {

public:

  shifted_lat() : GridCreator("^[Ss][Ll][Aa][Tt]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slat<gauss>" << "Shifted latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lat_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_lonlat_x : public GridCreator {

public:

  regular_lonlat_x() : GridCreator("^[Ll]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "L<nx>x<ny>" << "Regular longitude-latitude grid";
  }


  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::RegularLonLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} regular_lonlat_x_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lonlat_x : public GridCreator {

public:

  shifted_lonlat_x() : GridCreator("^[Ss]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "S<nx>x<ny>" << "Shifted longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLonLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lonlat_x_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lon_x : public GridCreator {

public:

  shifted_lon_x() : GridCreator("^[Ss][Ll][Oo][Nn]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slon<nx>x<ny>" << "Shifted longitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLon::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lon_x_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lat_x : public GridCreator {

public:

  shifted_lat_x() : GridCreator("^[Ss][Ll][Aa][Tt]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slat<nx>x<ny>" << "Shifted latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lat_x_;

}  // anonymous namespace

//---------------------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

