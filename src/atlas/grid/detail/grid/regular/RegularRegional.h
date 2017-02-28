#pragma once

#include "atlas/grid/detail/grid/Regular.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class RegularRegional: public Regular {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "regular_regional"; }

    RegularRegional( const Config& );

private:

    void setup(const util::Config& params);
    std::string shortName();

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
