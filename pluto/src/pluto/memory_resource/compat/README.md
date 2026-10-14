# pluto std::pmr memory_resource compatibility layer

Some C++ compilers and standard library combinations report C++17 support, but still do not provide the `<memory_resource>` header (or do not provide a usable implementation).

This directory provides a minimal API-compatible fallback for the parts of `std::pmr` that Pluto needs.

Important namespace note:

- Standard implementation: `std::pmr`
- Compatibility implementation in this directory: `pluto::compat`

Selection is handled in `pluto/memory_resource.h` via `PLUTO_HAVE_PMR`:

- When `PLUTO_HAVE_PMR = 1`, Pluto uses `std::pmr`.
- Otherwise, Pluto includes this compatibility implementation and uses `pluto::compat`.
- Currently by default `PLUTO_HAVE_PMR = 0`, and can be controlled via CMake configuration
  `-DENABLE_PMR=ON`
