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

    virtual std::string shortName() const { return "regular"; }
    virtual std::string gridType() const { return "regular"; }

    Regular(const util::Config& params);
    Regular();

    virtual eckit::Properties spec() const;

    size_t nlon() { return nlonmin(); }	// same for all latitudes
    double lon(size_t jlon) { return Structured::lon(0,jlon); } // same for all latitudes

protected:

    void setup();

protected:

    eckit::SharedPtr<spacing::Spacing> spacing_x_;
    eckit::SharedPtr<spacing::Spacing> spacing_y_;

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
