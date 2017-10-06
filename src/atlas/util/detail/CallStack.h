#pragma once

#include <list>
#include <iosfwd>
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace util {
namespace detail {

/// @class CallStack
/// Instances of CallStack can keep track of nested eckit::CodeLocations
class CallStack: public std::list< eckit::CodeLocation > {
private:

  using Base = std::list< eckit::CodeLocation >;

public:

  // Enable constructors from base class
  using Base::Base;

private:

  void print(std::ostream& out) const;

  friend std::ostream& operator<< ( std::ostream& out, const CallStack& p ) {
    p.print(out);
    return out;
  }

};

} // namespace detail
} // namespace util
} // namespace atlas
