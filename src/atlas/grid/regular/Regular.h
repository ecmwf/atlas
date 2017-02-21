#ifndef atlas_grid_regular_Regular_h
#define atlas_grid_regular_Regular_h

#include "atlas/grid/Structured.h"


namespace atlas {
namespace grid {
namespace regular {

class Regular: public Structured {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual eckit::Properties spec() const;

    size_t nlon()           { return nlonmin(); }               // same for all latitudes
    double lon(size_t jlon) { return Structured::lon(0,jlon); } // same for all latitudes

protected:

    Regular();

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
