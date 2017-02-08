#ifndef atlas_UniformSpacing_H
#define atlas_UniformSpacing_H

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

class UniformSpacing: public Spacing {


  public:

    // constructor
    UniformSpacing(const eckit::Parametrisation& p);

    // class name
    static std::string className() { return "atlas.UniformSpacing"; }
    static std::string spacing_type_str() {return "uniform";}

    void generate(size_t i, double &x) const;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas


#endif
