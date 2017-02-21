#ifndef atlas_grid_reduced_ReducedGaussian_h
#define atlas_grid_reduced_ReducedGaussian_h

#include "atlas/grid/Structured.h"

namespace atlas {
namespace grid {
namespace reduced {

class ReducedGaussian: public Structured {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "reduced_gaussian"; }

    ReducedGaussian(const util::Config& params);
    ReducedGaussian(const int N, const long pl[]);

    virtual eckit::Properties spec() const;

protected:

    ReducedGaussian();

    void setup(const size_t N, const long pl[]);

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
