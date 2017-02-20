#pragma once

#include <array>
#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

class CustomSpacing: public Spacing {

public:

    // constructor
    CustomSpacing(const eckit::Parametrisation& p);

    CustomSpacing(long N, const double x[], const std::array<double,2>& range = {-90.,90} );

    // class name
    static std::string className() { return "atlas.CustomSpacing"; }
    static std::string spacing_type_str() {return "custom";}

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
