
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <cstddef>

#include "pluto/pluto.h"

[[maybe_unused]] constexpr std::size_t kb = 1024;
[[maybe_unused]] constexpr std::size_t mb = 1024*kb;
[[maybe_unused]] constexpr std::size_t gb = 1024*mb;

void use_allocator() {
    // The allocator will use the default memory resource
    std::size_t size = 1*mb;
    pluto::allocator<std::byte> allocator;
    std::byte* array = allocator.allocate(size);
    allocator.deallocate(array, size);
}

void recursive_scope_push_pop(int recursion) {
    // Manually push new scope
    pluto::scope::push();

    pluto::TraceMemoryResource mr{"scope_"+std::to_string(recursion), pluto::new_delete_resource()};
    pluto::set_default_resource(&mr);
    use_allocator();
    if (recursion > 1) {
      recursive_scope_push_pop(recursion-1);
    }
    use_allocator();

    // Manually pop end scope
    pluto::scope::pop();
}

void recursive_scope_object(int recursion) {
    // Use scope object lifetime to push/pop scope
    pluto::scope scope;

    pluto::TraceMemoryResource mr{"scope_"+std::to_string(recursion), pluto::new_delete_resource()};
    pluto::set_default_resource(&mr);
    use_allocator();
    if (recursion > 1) {
      recursive_scope_object(recursion-1);
    }
    use_allocator();
}


int main(int argc, char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::TraceOptions::instance().enabled = true;

    // Just wrap default resource for tracing allocations
    pluto::TraceMemoryResource default_memory_resource{"default", pluto::get_default_resource()};
    pluto::set_default_resource(&default_memory_resource);

    use_allocator(); // should print default allocations

    recursive_scope_push_pop(3);

    use_allocator(); // should print default allocations

    recursive_scope_object(3);

    use_allocator(); // should print default allocations

    std::cout << "END" << std::endl;
}
