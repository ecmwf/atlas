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
static GridCreator::Registry *named_grids = 0;
static GridCreator::Registry *typed_grids = 0;


static pthread_once_t once = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    named_grids = new GridCreator::Registry();
    typed_grids = new GridCreator::Registry();
}

}  // anonymous namespace

//---------------------------------------------------------------------------------------------------------------------

const GridCreator::Registry& GridCreator::nameRegistry() {
  return *named_grids;
}

const GridCreator::Registry& GridCreator::typeRegistry() {
  return *typed_grids;
}


GridCreator::GridCreator( const std::string& type ) :
  names_(),
  type_(type) {
  pthread_once(&once, init);
  eckit::AutoLock<eckit::Mutex> lock(local_mutex);

  ASSERT(typed_grids->find(type_) == typed_grids->end());
  (*typed_grids)[type] = this;
}


GridCreator::GridCreator( const std::vector<std::string>& names ) :
  names_(names),
  type_() {
  pthread_once(&once, init);
  eckit::AutoLock<eckit::Mutex> lock(local_mutex);
  for( const std::string& name : names_ ) {
      ASSERT(named_grids->find(name) == named_grids->end());
      (*named_grids)[name] = this;
  }
}

GridCreator::GridCreator( const std::string& type, const std::vector<std::string>& names ) :
  names_(names),
  type_(type) {
  pthread_once(&once, init);
  eckit::AutoLock<eckit::Mutex> lock(local_mutex);

  for( const std::string& name : names_ ) {
      ASSERT(named_grids->find(name) == named_grids->end());
      (*named_grids)[name] = this;
  }

  ASSERT(typed_grids->find(type_) == typed_grids->end());
  (*typed_grids)[type] = this;
}

GridCreator::~GridCreator() {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    for( const std::string& name : names_ ) {
        ASSERT(named_grids->find(name) != named_grids->end());
        (*named_grids).erase(name);
    }

    if( not type_.empty() ) {
      ASSERT(typed_grids->find(type_) != typed_grids->end());
      (*typed_grids).erase(type_);
    }
}

const Grid::grid_t* GridCreator::create( const Grid::Config& config ) const {

  eckit::Factory<Grid::grid_t>& fact = eckit::Factory<Grid::grid_t>::instance();

  std::string name;
  if (config.get("name",name)) { // ignore any further configuration
      return create(name);
  }

  std::string type;
  if (config.get("type",type) && fact.exists(type)) {
      return fact.get(type).create(config);
  }

  if( name.size() ) {
    Log::error() << "name provided: " << name << std::endl;
  }
  if( type.size() ) {
    Log::error() << "type provided: " << type << std::endl;
  }
  if( name.empty() && type.empty() ) {
    throw eckit::BadParameter("no name or type in configuration",Here());
  } else {
    throw eckit::BadParameter("name or type in configuration don't exist",Here());
  }

  return nullptr;
}


bool GridCreator::match( const std::string& string, std::vector<std::string>& matches, int &id ) const {
  id = 0;
  for( const std::string& name : names_ ) {
    if( Regex(name).match(string,matches) )
      return true;
    ++id;
  }
  return false;
}

std::string GridCreator::type() const { return type_; }

std::ostream& operator<<( std::ostream& os, const GridCreator& g ) {
  g.print(os);
  return os;
}

//---------------------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

