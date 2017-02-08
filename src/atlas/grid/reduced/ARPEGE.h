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

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "arpege"; }


    ARPEGE(const util::Config& params);

protected:

    void setup(size_t N);

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
