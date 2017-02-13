/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/grids.h"

#include <regex.h>
#include "eckit/parser/Tokenizer.h"
#include "eckit/utils/Translator.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"


namespace atlas {
namespace grid {


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
        printf("This regular expression didn't compile: \"%s\"", regex.c_str());

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
    bool match(const std::string& string) {
        std::vector<std::string> substr;
        return regex_match_impl(string,regex_,substr,false,use_case_);
    }
    bool match(const std::string& string, std::vector<std::string>& substr) {
        return regex_match_impl(string,regex_,substr,true,use_case_);
    }
  private:
    std::string regex_;
    bool use_case_;
};

//=====================================================

Grid* grid_from_uid(const std::string& uid) {
    if (eckit::Factory<Grid>::instance().exists(uid)) {
        return Grid::create(util::Config("grid_type", uid));
    } else {
        Regex classical_gaussian  ("^[Nn]([0-9]+)$");
        Regex octahedral_gaussian ("^[Oo]([0-9]+)$");
        Regex regular_gaussian    ("^[Ff]([0-9]+)$");
        Regex regular_lonlat      ("^[Ll]([0-9]+)$");
        Regex shifted_lonlat      ("^[Ss]([0-9]+)$");
        Regex shifted_lon         ("^[Ss][Ll][Oo][Nn]([0-9]+)$");
        Regex shifted_lat         ("^[Ss][Ll][Aa][Tt]([0-9]+)$");
        Regex regular_lonlat_x    ("^[Ll]([0-9]+)x([0-9]+)$");
        Regex shifted_lonlat_x    ("^[Ss]([0-9]+)x([0-9]+)$");
        Regex shifted_lon_x       ("^[Ss][Ll][Oo][Nn]([0-9]+)x([0-9]+)$");
        Regex shifted_lat_x       ("^[Ss][Ll][Aa][Tt]([0-9]+)x([0-9]+)$");

        util::Config gridparams;
        eckit::Translator<std::string,int> to_int;
        std::vector<std::string> matches;
        if( classical_gaussian.match(uid,matches) ) {
            try {
                int N = to_int(matches[0]);
                gridparams.set("grid_type", reduced::ClassicGaussian::grid_type_str());
                gridparams.set("N",N);
                return Grid::create( gridparams );
            } catch( const eckit::BadParameter& e ) {
                throw eckit::BadParameter("Grid ["+uid+"] does not exist.\n"+e.what(),Here());
            }
            return 0;
        } else if( octahedral_gaussian.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", reduced::OctahedralGaussian::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( regular_gaussian.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", regular::RegularGaussian::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( regular_lonlat.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", regular::RegularLonLat::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( shifted_lonlat.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", regular::ShiftedLonLat::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( shifted_lon.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", regular::ShiftedLon::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( shifted_lat.match(uid,matches) ) {
            int N = to_int(matches[0]);
            gridparams.set("grid_type", regular::ShiftedLat::grid_type_str());
            gridparams.set("N",N);
            return Grid::create( gridparams );
        } else if( regular_lonlat_x.match(uid,matches) ) {
            int nlon = to_int(matches[0]);
            int nlat = to_int(matches[1]);
            gridparams.set("grid_type", regular::RegularLonLat::grid_type_str());
            gridparams.set("nlon",nlon);
            gridparams.set("nlat",nlat);
            return Grid::create( gridparams );
        } else if( shifted_lonlat_x.match(uid,matches) ) {
            int nlon = to_int(matches[0]);
            int nlat = to_int(matches[1]);
            gridparams.set("grid_type", regular::ShiftedLonLat::grid_type_str());
            gridparams.set("nlon",nlon);
            gridparams.set("nlat",nlat);
            return Grid::create( gridparams );
        } else if( shifted_lon_x.match(uid,matches) ) {
            int nlon = to_int(matches[0]);
            int nlat = to_int(matches[1]);
            gridparams.set("grid_type", regular::ShiftedLon::grid_type_str());
            gridparams.set("nlon",nlon);
            gridparams.set("nlat",nlat);
            return Grid::create( gridparams );
        } else if( shifted_lat_x.match(uid,matches) ) {
            int nlon = to_int(matches[0]);
            int nlat = to_int(matches[1]);
            gridparams.set("grid_type", regular::ShiftedLat::grid_type_str());
            gridparams.set("nlon",nlon);
            gridparams.set("nlat",nlat);
            return Grid::create( gridparams );
        } else {
            eckit::Tokenizer tokenize(".");
            std::vector<std::string> tokens;
            tokenize(uid,tokens);
            std::string grid_type = tokens[0];
            if( tokens.size() > 1) {
                gridparams.set("grid_type",grid_type);
                if( tokens[1][0] == 'N' ) {
                    std::string Nstr(tokens[1],1,tokens[1].size()-1);
                    int N = to_int(Nstr);
                    gridparams.set("N",N);
                } else {
                    std::vector<std::string> lonlat;
                    eckit::Tokenizer tokenize_lonlat("x");
                    tokenize_lonlat(tokens[1],lonlat);
                    if( lonlat.size() > 1 ) {
                        int nlon = to_int(lonlat[0]);
                        int nlat = to_int(lonlat[1]);
                        gridparams.set("nlon",nlon);
                        gridparams.set("nlat",nlat);
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
void load_grid() {
    eckit::ConcreteBuilderT1<Grid,CONCRETE> builder("tmp");
}

void load() {
    Log::debug(2) << "Loading library [atlas::grid]" << std::endl;

    // We have to touch all classes we want to register for static linking.

    load_grid<Unstructured>();
    load_grid<CustomStructured>();
    load_grid<reduced::ReducedGaussian>();
    load_grid<regular::RegularGaussian>();
    load_grid<reduced::ClassicGaussian>();
    load_grid<reduced::OctahedralGaussian>();
    load_grid<reduced::ARPEGE>();
    load_grid<reduced::ReducedLonLat>();
    load_grid<regular::RegularLonLat>();
    load_grid<regular::ShiftedLonLat>();
    load_grid<regular::ShiftedLon>();
    load_grid<regular::ShiftedLat>();

}


extern "C"
{

    void atlas__grids__load() {
        atlas::grid::load();
    }

}


} // namespace grid
} // namespace atlas

