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

    RegularRegional(const util::Config& params);
    RegularRegional();

    eckit::Properties spec() const;
    
  protected:

    void setup(const util::Config& params);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
