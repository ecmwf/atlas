#ifndef atlas_UniformSpacing_H
#define atlas_UniformSpacing_H

#include <array>
#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

class UniformSpacing: public Spacing {

public:

    // constructor
    UniformSpacing( const eckit::Parametrisation& p );

    UniformSpacing( double centre, double step, long N, bool endpoint=true );

    UniformSpacing( std::array<double,2> min_max, long N, bool endpoint=true );

    // class name
    static std::string className() { return "atlas.UniformSpacing"; }
    static std::string spacing_type_str() {return "uniform";}

protected:

    // points are equally spaced between xmin and xmax
    // Depending on value of endpoint, the spacing will be different
    void setup(double min, double max, long N, bool endpoint);

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas


#endif
