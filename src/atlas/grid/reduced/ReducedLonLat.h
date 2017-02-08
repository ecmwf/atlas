#ifndef atlas_grid_reduced_ReducedLonLat_h
#define atlas_grid_reduced_ReducedLonLat_h

#include "atlas/grid/Structured.h"

namespace atlas {
namespace grid {
namespace reduced {

class ReducedLonLat: public Structured {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "reduced_lonlat"; }

    ReducedLonLat(): Structured() {}
    ReducedLonLat(const util::Config& params);

    eckit::Properties spec() const;

protected:

    void setup(size_t nlat, long pl[]);

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
