#ifndef atlas_grid_regular_RegularRegional_h
#define atlas_grid_regular_RegularRegional_h

#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class RegularRegional: public Regular {

  public:

    static std::string grid_type_str();

    static std::string className();
    
    virtual std::string shortName() const;

    RegularRegional(const util::Config& params);
    RegularRegional();
    
  protected:

    void setup(const util::Config& params);
    std::string shortName();

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
