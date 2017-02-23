#ifndef atlas_Spacing_H
#define atlas_Spacing_H

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

    typedef const eckit::Parametrisation& ARG1;
    typedef eckit::BuilderT1<Spacing> builder_t;
    typedef std::vector<double>::const_iterator const_iterator;
    using Interval = std::array<double,2>;

public:

    static Spacing* create(const eckit::Parametrisation& params);

    static std::string className() {return "atlas.Spacing";}
    static std::string spacing_type_str() {return "spacing";}

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

protected:

    std::vector<double> x_;
    double min_;
    double max_;

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas


#endif
