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

    static std::string static_type();

    virtual std::string name() const;

    RegularRegional( const Config& );

private:

    void setup(const util::Config& params);
    std::string name();

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
