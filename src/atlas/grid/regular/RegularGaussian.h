#ifndef atlas_grid_regular_RegularGaussian_h
#define atlas_grid_regular_RegularGaussian_h

#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class RegularGaussian: public Regular {

  public:

    static std::string grid_type_str();

    static std::string className();
    
    virtual std::string shortName() const;

    RegularGaussian(const util::Config& params);

  protected:

    void setup();

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
