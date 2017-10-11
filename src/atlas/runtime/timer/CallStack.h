#pragma once

#include <list>

namespace eckit { class CodeLocation; }

namespace atlas {
namespace runtime {
namespace timer {

/// @class CallStack
/// Instances of CallStack can keep track of nested eckit::CodeLocations
class CallStack {

public:

  using const_iterator = std::list<size_t>::const_iterator;
  using const_reverse_iterator = std::list<size_t>::const_reverse_iterator;

public:

  void push_front( const eckit::CodeLocation& );
  void pop_front();

  const_iterator begin() const { return stack_.begin(); }
  const_iterator end()   const { return stack_.end(); }

  const_reverse_iterator rbegin() const { return stack_.rbegin(); }
  const_reverse_iterator rend()   const { return stack_.rend(); }

  size_t hash() const;
  size_t size() const { return stack_.size(); }

private:

  std::list<size_t> stack_;
  mutable size_t hash_{0};

};

} // namespace timer
} // namespace runtime
} // namespace atlas
