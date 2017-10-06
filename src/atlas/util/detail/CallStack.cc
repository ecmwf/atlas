
#include <iostream>

#include "CallStack.h"

namespace atlas {
namespace util {
namespace detail {

void CallStack::print(std::ostream& out) const {
    const_iterator it = begin();
    out << "\nstack:\n";
    for( size_t i=0; i<size(); ++i, ++it ) {
      out << i << ": " << it->file() << " +"<< it->line() << '\n';
    }
    out << "\n";
    out << std::flush;
}

} // namespace detail
} // namespace util
} // namespace atlas
