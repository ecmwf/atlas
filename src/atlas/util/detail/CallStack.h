#pragma once

#include <list>
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace util {
namespace detail {

/// @class CallStack
/// Instances of CallStack can keep track of nested eckit::CodeLocations
class CallStack {

public:
  using const_iterator = std::list<size_t>::const_iterator;
  using const_reverse_iterator = std::list<size_t>::const_reverse_iterator;

public:

  void push_front( const eckit::CodeLocation& );
  void pop_front();

  const eckit::CodeLocation& loc() const { return loc_; }

  const_iterator begin() const { return stack_.begin(); }
  const_iterator end()   const { return stack_.end(); }

  const_reverse_iterator rbegin() const { return stack_.rbegin(); }
  const_reverse_iterator rend()   const { return stack_.rend(); }

  size_t hash() const;
  size_t size() const { return stack_.size(); }

private:

  eckit::CodeLocation loc_;
  std::list<size_t> stack_;
  mutable size_t hash_{0};

};

} // namespace detail
} // namespace util
} // namespace atlas
