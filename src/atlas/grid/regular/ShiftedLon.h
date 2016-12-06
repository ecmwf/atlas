#ifndef atlas_grid_regular_ShiftedLon_h
#define atlas_grid_regular_ShiftedLon_h

#include "atlas/grid/regular/RegularLonLat.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class ShiftedLon: public RegularLonLat {

  public:

    static std::string grid_type_str();

    static std::string className();

    ShiftedLon(const util::Config& params);

    eckit::Properties spec() const;
    
  protected:

    void setup();

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
