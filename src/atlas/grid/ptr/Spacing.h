#pragma once

#include "atlas/grid/spacing/Spacing.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {
namespace ptr {

//---------------------------------------------------------------------------------------------------------------------

class Spacing {

public:

  using const_iterator = atlas::grid::spacing::Spacing::const_iterator;
  using Interval = atlas::grid::spacing::Spacing::Interval;

public:

  Spacing();
  Spacing( atlas::grid::spacing::Spacing* );
  Spacing( const eckit::Parametrisation& );

  operator atlas::grid::spacing::Spacing*() { return spacing_; }


  size_t size() const { return spacing_->size(); }

  double operator[](size_t i) const { return spacing_->operator[](i); }

  const_iterator begin() const { return spacing_->begin(); }
  const_iterator end()   const { return spacing_->end();   }

  const double& front() const { return spacing_->front(); }
  const double& back()  const { return spacing_->back();  }

  Interval interval() const { return spacing_->interval(); }

  const double min() const { return spacing_->min(); }
  const double max() const { return spacing_->max(); }


private:

  atlas::grid::spacing::Spacing* spacing_;
};

//---------------------------------------------------------------------------------------------------------------------

}
}
}
