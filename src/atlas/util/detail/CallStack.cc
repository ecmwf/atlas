
#include <functional>

#include "CallStack.h"

namespace atlas {
namespace util {
namespace detail {

void CallStack::push_front( const eckit::CodeLocation& loc ) {
  stack_.push_front( std::hash<std::string>{}(loc.asString()) );
  loc_ = loc;
}

void CallStack::pop_front() {
  stack_.pop_front();
}

size_t CallStack::hash() const {
  if( hash_ ) return hash_;
  for( auto h : stack_ ) {
    hash_ ^= (h << 1);
  }
  return hash_;
}

} // namespace detail
} // namespace util
} // namespace atlas
