#ifndef atlas_grid_regular_GlobalLonLat_h
#define atlas_grid_regular_GlobalLonLat_h

#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class GlobalLonLat: public Regular {

  public:

    static std::string grid_type_str();

    static std::string className();
		
		virtual std::string shortName() const;
		virtual std::string gridType() const { return "global_lonlat"; }
		
    GlobalLonLat();
    GlobalLonLat(const util::Config& params);
    GlobalLonLat(long nlon, long nlat);
    GlobalLonLat(long N);

    eckit::Properties spec() const;
    
    bool isShiftedLon() const { return shiftLon_; }
    bool isShiftedLat() const { return shiftLat_; }
    
  protected:

    void setup(long nlon, long nlat);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);
    
    bool shiftLon_;
    bool shiftLat_;

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
