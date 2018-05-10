
#include "CallStack.h"

#include <functional>

#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace runtime {
namespace trace {

void CallStack::push_front( const eckit::CodeLocation& loc, const std::string& id ) {
    stack_.push_front( std::hash<std::string>{}( loc.asString() + id ) );
}

void CallStack::pop_front() {
    stack_.pop_front();
}

size_t CallStack::hash() const {
    if ( hash_ ) return hash_;
    for ( auto h : stack_ ) {
        hash_ ^= ( h << 1 );
    }
    return hash_;
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
