#ifndef atlas_grid_reduced_ARPEGE_h
#define atlas_grid_reduced_ARPEGE_h

#include "atlas/grid/reduced/ClassicGaussian.h"

namespace atlas {
namespace grid {
namespace reduced {

class ARPEGE: public ClassicGaussian {

  public:

    static std::string grid_type_str();

    static std::string className();

    ARPEGE(const util::Config& params);

    eckit::Properties spec() const;
    
    virtual const domain::Domain& domain() const { return *domain_; }
    
  protected:

    void setup(size_t N);

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
