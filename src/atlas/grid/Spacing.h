#pragma once

#include "atlas/grid/detail/spacing/Spacing.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Spacing {

public:

  using const_iterator = atlas::grid::spacing::Spacing::const_iterator;
  using Interval = atlas::grid::spacing::Spacing::Interval;

public:

  Spacing();
  Spacing( const Spacing& );
  Spacing( const atlas::grid::spacing::Spacing* );
  Spacing( const eckit::Parametrisation& );

  operator bool() const { return spacing_; }

  operator const atlas::grid::spacing::Spacing*() { return spacing_.get(); }

  size_t size() const { return spacing_.get()->size(); }

  double operator[](size_t i) const { return spacing_.get()->operator[](i); }

  const_iterator begin() const { return spacing_.get()->begin(); }
  const_iterator end()   const { return spacing_.get()->end();   }

  const double& front() const { return spacing_.get()->front(); }
  const double& back()  const { return spacing_.get()->back();  }

  Interval interval() const { return spacing_.get()->interval(); }

  const double min() const { return spacing_.get()->min(); }
  const double max() const { return spacing_.get()->max(); }

  std::string type() const { return spacing_.get()->type(); }


private:

  eckit::SharedPtr<const atlas::grid::spacing::Spacing> spacing_;
};

//---------------------------------------------------------------------------------------------------------------------

}
}
