#pragma once

#include <vector>
#include <array>
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"
#include "atlas/util/Config.h"

namespace eckit {
  class Parametrisation;
}

namespace atlas {
namespace grid {
namespace spacing {

class Spacing : public eckit::Owned {

public:

    using const_iterator = std::vector<double>::const_iterator;
    using Interval       = std::array<double,2>;

    using ARG1      = const eckit::Parametrisation&;
    using builder_t = eckit::BuilderT1<Spacing>;
    static std::string className() {return "atlas.Spacing";}

public:

    static Spacing* create(const eckit::Parametrisation& params);

    virtual std::string type() const =0;

    double operator[](size_t i) const { return x_[i]; }

    size_t size() const { return x_.size(); }

    const double* data() const { return x_.data(); }

    const_iterator begin() const { return x_.begin(); }
    const_iterator end()   const { return x_.end();   }

    const double& front() const { return x_.front(); }
    const double& back()  const { return x_.back();  }

    Interval interval() const { return {min_,max_}; }

    const double min() const { return min_; }
    const double max() const { return max_; }
    
    virtual eckit::Properties spec() const =0;

protected:

    std::vector<double> x_;
    double min_;
    double max_;

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
