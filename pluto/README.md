Pluto
=====

In Greek mythology, Plouton is the god of the underworld,
often associated with wealth and the riches found underground,
such as precious metals and minerals. 
This name was later Latinized to Pluto in Roman mythology.

What is it?
===========

Pluto contains high-level abstractions for memory resource management and offloading data.
The memory resource management is based on and compatible with C++17 `std::pmr::memory_resource` and `std::pmr::polymorphic_allocator`,
and is extended with asynchronous (de)allocation methods.

GPU specific memory resources are available, delegating to the low-level hic (HIP/CUDA abstraction) library to support GPUs.

Pluto can be used and configured both from C++ and Fortran.

The concepts
============

### pluto::memory_resource

The `pluto::memory_resource` abstract class is an alias for
[`std::pmr::memory_resource`](https://en.cppreference.com/w/cpp/memory/memory_resource)
and provides following noteworthy member functions:
```c++
void* allocate(std::size_t bytes, std::size_t alignment) {
    return do_allocate(bytes, alignment);
}
void deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    do_deallocate(ptr, bytes, alignment);
}
```
Concrete implementations deriving from `pluto::memory_resource` must implement these functions:
```c++
void* do_allocate(std::size_t bytes, std::size_t alignment) override;
void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
```

Pluto provides 8 predefined concrete implementations:
| memory_resource         | memory_pool_resource         |
|-------------------------|------------------------------|
| pluto::host_resource    | pluto::host_pool_resource    |
| pluto::device_resource  | pluto::device_pool_resource  |
| pluto::pinned_resource  | pluto::pinned_pool_resource  |
| pluto::managed_resource | pluto::managed_pool_resource |

These predefined memory resources have unlimited lifetime and have memory tracking and tracing capability.  \
See [Predefined pluto::memory\_resources](#predefined-plutomemory_resources) below for details on each memory resource.

For convenience, pluto also provides aliases to two predefined standard library `std::pmr::memory_resources`:
- pluto::new_delete_resource -> [std::pmr::new_delete_resource](https://en.cppreference.com/w/cpp/memory/new_delete_resource)
- pluto::null_memory_resource -> [std::pmr::null_memory_resource](https://en.cppreference.com/w/cpp/memory/null_memory_resource)

#### Example:
```C++
double* data;
std::size_t bytes     = 10 * sizeof(double);
std::size_t alignment = 64;
pluto::memory_resource* mr = pluto::host_resource();
double* data = (double*) mr->allocate(bytes, alignment);
mr->deallocate(data, bytes, alignment);
```

### pluto::async_memory_resource

The `pluto::async_memory_resource` extends `pluto::memory_resource` with asynchronous allocation and deallocation features. The asynchronous argument is a `pluto::stream_view`, which implements a `cudaStream` or `hipStream`. \
The extra member functions are:
```c++
void* allocate_async(std::size_t bytes, std::size_t alignment, pluto::stream_view stream) {
    return do_allocate_async(bytes, alignment, stream);
}
void deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, pluto::stream_view stream) {
    do_deallocate_async(ptr, bytes, alignment, stream);
}
```
Concrete implementations deriving from `pluto::async_memory_resource` then further implement these functions:
```c++
void* do_allocate_async(std::size_t bytes, std::size_t alignment, pluto::stream_view stream) override;
void do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, pluto::stream_view stream) override;
```

### pluto::allocator

The `pluto::allocator<T>` extends [`std::pmr::polymorphic_allocater`](https://en.cppreference.com/w/cpp/memory/polymorphic_allocator)
which implements all functions required of a C++ [Allocator](https://en.cppreference.com/w/cpp/named_req/Allocator) to be given to [AllocatorAwareContainers](https://en.cppreference.com/w/cpp/named_req/AllocatorAwareContainer) such as e.g. `std::vector`, `std::map`, `std::set`, `std::list`, `std::string`.
It internally uses a `pluto::memory_resource*` for allocation and deallocation.
The noteworthy functions are:
```c++
/// Constructor without arguments; a configurable default pluto::memory_resource will be used
/// This default can be set with `std::pmr::set_default_resource()` or `pluto::set_default_resource()`
pluto::allocator<T>();

/// Constructor using a given memory_resource. Note this is compatible with any third-party `std::pmr::memory_resource`
pluto::allocator<T>(pluto::memory_resource*);

/// Return an new allocated array with `size` number of elements. (Not bytes unlike pluto::memory_resource)
T* allocate(std::size_t size);

/// Deallocate a given array with `size` number of elements. (Not bytes unlike pluto::memory_resource)
void deallocate(T* ptr, std::size_t size);

/// Return an new allocated array with `size` number of elements. (Not bytes unlike pluto::memory_resource)
T* allocate_async(std::size_t size, pluto::stream_view stream);

/// Deallocate a given array with `size` number of elements. (Not bytes unlike pluto::memory_resource)
void deallocate_async(T* ptr, std::size_t size, pluto::stream_view stream);
```

When `(de)allocate_async` is used with a `memory_resource` that does not derive from `pluto::async_memory_resource`, then `(de)allocate` will be used instead.
The functions `(de)allocate`, `(de)allocate_async` also have overrides with a first argument `std::string_view label`, that can be used for tracing memory (de)allocations if the used concrete memory resources supports it.

#### Examples:

- Use allocator to allocate array of 10 elements using pre-defined `pluto::host_resource()` \
  C++:
  ```c++
  pluto::allocator<double> alloc(pluto::host_resource());
  double* data = alloc.allocate(10);
  alloc.deallocate(data, 10);
  ```
  Fortran:
  ```fortran
  real(8), pointer :: array1d(:)
  type(pluto_allocator) :: alloc
  alloc = pluto%make_allocator(pluto%host_resource())
  call alloc%allocate(array1d, [10])
  call alloc%deallocate(array1d)
  ```
- Use memory pool in allocator-aware type using pre-defined `pluto::host_pool_resource()`
  ```c++
  std::vector<double, pluto::allocator<double>> vector(pluto::host_pool_resource());
  vector.resize(10);
  ```
- The latter can be done via the `std::pmr::vector` as well due to the `std::pmr::memory_resource` compatibility:
   ```c++
   std::pmr::vector<double> vector(pluto::host_pool_resource());
   vector.resize(10);
   ```
- We don't need to explicitely add the `memory_resource` in the `std::pmr::vector` constructor,
  when setting the default beforehand:
   ```c++
   std::pmr::set_default_resource(pluto::host_pool_resource());
   std::pmr::vector<double> vector;
   vector.resize(10);
   ```

### pluto::{host,device} namespace

In namespaces `pluto::host` and `pluto::device`, pluto manages defaults per memory space, independently from `std::pmr::{get,set}_default_resource()`.

Following functions exist for C++
```c++
pluto::host::set_default_resource(pluto::memory_resource*);
pluto::host::get_default_resource() -> pluto::memory_resource*;
pluto::device::set_default_resource(pluto::memory_resource*);
pluto::device::get_default_resource() -> pluto::memory_resource*;
```
Following routines exist for Fortran (pseudocode):
```
pluto%host%set_default_resource( type(pluto_memory_resource) )
pluto%host%get_default_resource() -> type(pluto_memory_resource)
pluto%device%set_default_resource( type(pluto_memory_resource) )
pluto%device%get_default_resource() -> type(pluto_memory_resource)
```

Following C++ classes exist that extend `pluto::allocator<T>`:
```c++
pluto::host::allocator<T>
pluto::device::allocator<T>
```
The only difference with `pluto::allocator` is that the default constructor won't use `std::pmr::get_default_resource()`, \
but rather `pluto::{host,device}::get_default_resource()`.

In Fortran you would create allocators that use the memory space specific defaults via:
```fortran
type(pluto_allocator) :: host_alloc, device_alloc
host_alloc   = pluto%host%make_allocator()
device_alloc = pluto%device%make_allocator()
```

The initial value returned by `pluto::{host,device}::get_default_resource()` is respectively `pluto::host_resource()` and `pluto::device_resource()` unless specified otherwise via environment variables:
```sh
export PLUTO_HOST_MEMORY_RESOURCE=pluto::pinned_pool_resource
export PLUTO_DEVICE_MEMORY_RESOURCE=pluto::device_pool_resource
```

Predefined pluto::memory_resources
----------------------------------

Pluto provides a number of predefined concrete `pluto::memory_resources` via accessor
 returning `pluto::memory_resource*` in C++ or `type(pluto_memory_resource)` in Fortran.
They have C++ and Fortran accessor functions and are as well registered by name:
- **pluto::new_delete_resource** \
   Alias to `std::pmr::new_delete_resource`, using C++ new and delete.

- **pluto::null_memory_resource()** \
   Alias to `std::pmr::null_memory_resource`, throwing exception when used.

- **pluto::host_resource** \
   Allocates host CPU memory aligned to 256 bytes.

- **pluto::host_pool_resource** \
   A memory pool based on pluto::host_resource

- **pluto::pinned_resource** \
   Allocates host-pinned (a.k.a. page-locked) CPU memory aligned to 256 bytes.

- **pluto::pinned_pool_resource** \
   A memory pool based on pluto::pinned_resource

- **pluto::device_resource** \
   A `pluto::async_memory_resource` that allocates device resident memory. \
   Internally this uses `cudaMalloc` or `hipMalloc` for allocate and
   `cudaMallocAsync` or `hipMallocAsync` for allocate_async

- **pluto::device_pool_resource** \
   A memory pool based on pluto::device_resource

- **pluto::managed_resource** \
   Allocates UVM a.k.a. managed memory accessible from both host and device.\
   Internally this uses `cudaMallocManaged` or `hipMallocManaged`

- **pluto::managed_pool_resource** \
   A memory pool based on pluto::managed_resource

## Data transfer

### pluto::memcpy\_{host,device}\_to\_{device,host}

These functions work with void* pointers and bytes arguments. \
An optional `pluto::stream_view` provides async data transfers

### pluto::copy\_{host,device}\_to\_{device,host}

These functions work with templated T* pointers, and size (number of elements) arguments. \
An optional `pluto::stream_view` provides async data transfers


### memcpy\_{host,device}\_to\_{device,host}\_2D and  copy\_{host,device}\_to\_{device,host}\_2D

Like the above but for discontiguous slices, this is useful for the atlas::MultiField


## Tracing and tracking memory

The pluto predefined memory resources have tracking and tracing capability.
To enable tracing, e.g. for debugging, set environment variable `PLUTO_TRACE=1`.
Fine control is also possible programaticaly:
```c++
bool previous_status = pluto::trace::enable(true);
// ... do stuff ...
pluto::trace::enable(previous_status);
```
The trace output gets written to `pluto::trace::out` stream, which defaults to `std::cout`. This can be modified, e.g.
```c++
std::stringstream pluto_trace_stream;
pluto::trace::set(pluto_trace_stream);
```

A memory usage report can be obtained, which reports on the use of each of the pluto predefined memory resources.
Other user-defined memory resources are not taken into consideration.
```
pluto::memory::report() -> std::string
```

# Real Examples

See examples subdirectory on how to use Pluto.
