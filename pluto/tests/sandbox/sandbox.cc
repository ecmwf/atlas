
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>
//#include <span>
#include <iostream>
#include <memory>
#include <new>

#include "pluto/pluto.h"

#include "kernel.h"

using namespace pluto;
using namespace pluto;

struct CustomType {
  HIC_HOST_DEVICE
  CustomType(bool verbose = false) :
    verbose_(verbose) {
    a = 999;
    b = 4;
    if(verbose_) {
#if HIC_DEVICE_COMPILE
    printf("    + CustomType(DEVICE) a=%d, b=%d\n", a, b);
#else
    printf("    + CustomType(HOST) a=%d, b=%d\n", a, b);
#endif
    }
  }
  HIC_HOST_DEVICE
  ~CustomType(){
    if( verbose_ ) {
#if HIC_DEVICE_COMPILE
    printf("    - ~CustomType(DEVICE) a=%d, b=%d\n", a, b);
#else
    printf("    - ~CustomType(HOST) a=%d, b=%d\n", a, b);
#endif
    }
  }
  HIC_HOST_DEVICE
  void compute() {
     a = 10;
     b = 11;
  }
  int a;
  int b;
  bool verbose_;
};

namespace pluto::host {
template<typename T>
class Array {
public:
    using value_type = T;

    Array() :
      resource_(get_default_resource()) {
    }

    Array(std::size_t size, std::size_t alignment=alignof(value_type)) : Array() {
      alignment_ = alignment;
      allocate(size);
    }
    ~Array() {
      deallocate();
    }

    //operator std::span<value_type,std::dynamic_extent>() { return std::span<value_type,std::dynamic_extent>(ptr_,size_); }

    void assign(std::size_t size, T value) {
      if(size_ == 0) {
        allocate(size);
      }
      //std::span<T> span(ptr_,size_);
      auto* span = ptr_; 
      for(size_t j=0; j<size; ++j) {
        span[j] = value;
      }
    }

    void allocate(std::size_t size) {
      size_ = size;
      if (size_ > 0) {
        //ptr_ = alloc_.allocate(size_);
        ptr_ = (value_type*)resource_->allocate(size_*sizeof(value_type), alignment_);
      }
      if (! is_aligned( ptr_, alignment_ )) {
        std::cout << "assert(is_aligned(ptr_, alignment_)) failed" << std::endl;
      }
    }

    void deallocate() {
      if (size_) {
        //alloc_.deallocate(ptr_, size_);
        resource_->deallocate(ptr_, size_*sizeof(value_type), alignment_);
      }
    }

private:
    memory_resource* resource_;
    //std::pmr::polymorphic_allocator<T> alloc_;
    value_type* ptr_{nullptr};
    std::size_t size_{0};
    std::size_t alignment_{256};
    //std::size_t alignment_{alignof(value_type)};
};
}

namespace pluto::device {
template<typename T>
class Array {
public:
    using value_type = T;

    Array() :
      Array(0, pluto::get_default_stream()) {}

    Array(std::size_t size) : 
      Array(size, pluto::get_default_stream()) {}

    Array(const pluto::Stream& stream) :
      Array(0, stream) {}

    Array(std::size_t size, const pluto::Stream& stream) : 
      resource_(pluto::device::get_default_resource()),
      stream_(stream),
      alloc_(resource_) {
      allocate(size);
    }

    ~Array() {
      deallocate();
    }

    //operator std::span<value_type>() { return std::span<value_type>(ptr_,size_); }

    void assign(std::size_t size, T value) {
      if(size_ == 0) {
        allocate(size);
      }
      //std::span<T> span(ptr_,size_);
      auto* span = ptr_;
      for(size_t j=0; j<size; ++j) {
        span[j] = value;
      }
    }

    void allocate(std::size_t size) {
      size_ = size;
      if (size_ > 0) {
          alloc_.allocate_async(size, stream_);
          // scoped_default_stream default_stream{stream_};
          // ptr_ = (value_type*)resource_->allocate(size_*sizeof(value_type),alignment_);
          if (! is_aligned( ptr_, alignment_ )) {
            std::cout << "assert(is_aligned(ptr_, alignment_)) failed" << std::endl;
          }
      }
    }

    void deallocate() {
      if (size_) {
          alloc_.deallocate_async(ptr_, size_, stream_);
          // scoped_default_stream default_stream{stream_};
          // resource_->deallocate(ptr_, size_*sizeof(value_type),alignment_);
        }
    }

private:
    pluto::allocator<value_type> alloc_;
    memory_resource* resource_;
    const pluto::Stream& stream_;
    value_type* ptr_{nullptr};
    std::size_t size_{0};
    static constexpr std::size_t alignment_ = pluto::default_alignment();
};
}

#if 0
template <typename T>
std::ostream& operator<< (std::ostream &os, std::span<T> v)
{
    os << "[ "; // <--
    for (const auto &i : v) {
        os << i << " ";
    }
    os << "]";
    return os;
}

template <typename T>
void print_vec(std::span<T> vec) {
  std::cout << vec << std::endl;
}
#endif

auto bytes_offset = [] (void* buffer, void* ptr) {
  return std::distance(reinterpret_cast<std::byte*>(buffer), reinterpret_cast<std::byte*>(ptr));
};

  template< typename T >
  // using vector = std::pmr::vector<T>;
  using vector = host::Array<T>;


int main(int argc, char* argv[]) {

  // host::register_resource("heap",
  //   std::make_unique<TraceMemoryResource>(std::cout, "heap",std::pmr::new_delete_resource())
  // );
  if (argc > 1) {
    TraceOptions::instance().enabled = std::atoi(argv[1]);
  }
#if 0
  std::size_t alignment = 64;
  std::size_t bytes = 32;
  // std::pmr::memory_resource* mr = std::pmr::new_delete_resource();
  auto new_delete_resource = std::make_unique<TraceMemoryResource>("new_delete_resource", std::pmr::new_delete_resource());

  auto pinned_new_delete =  std::make_unique<TraceMemoryResource>("pinned_new_delete", std::make_unique<PinnedMemoryResource>(new_delete_resource.get()));

  auto mr = std::make_unique<TraceMemoryResource>("pool", std::make_unique<host::AlignedMemoryPoolResource>(pinned_new_delete.get()));

  std::pmr::set_default_resource(new_delete_resource.get());

  if (false ){
    //auto mr = std::pmr::get_default_resource();
    void* data;
    
    data = mr->allocate(bytes, alignment);
    mr->deallocate(data, bytes, alignment);

    data = mr->allocate(bytes, alignment);
    mr->deallocate(data, bytes, alignment);

    data = mr->allocate(bytes*100, alignment);
    mr->deallocate(data, bytes*100, alignment);

    data = mr->allocate(bytes*200, alignment);
    mr->deallocate(data, bytes*200, alignment);

    type(atlas_MemoryResource) :: alloc
    alloc = get_registered_resource("heap")

    call kernel( .. , alloc )
        real(8), pointer :: tmp(:,:,:)

        if( not present(alloc))
            alloc = get_default_alloc()
        endif
        call alloc%allocate(tmp,shape)

  }

//  auto dmr = pluto::managed_resource();
  auto dmr = std::make_unique<TraceMemoryResource>("device_pool",std::make_unique<host::AlignedMemoryPoolResource>(pluto::device_resource()));

  {
    void* data;

    //pluto:::stream s{1};

    data = dmr->allocate(bytes, alignment);
    dmr->deallocate(data, bytes, alignment);

    data = dmr->allocate(bytes, alignment);
    dmr->deallocate(data, bytes, alignment);

    //dmr->release();

    data = dmr->allocate(bytes, alignment);
    dmr->deallocate(data, bytes, alignment);

    data = dmr->allocate(bytes, alignment);
    dmr->deallocate(data, bytes, alignment);

  }

  std::cout << "\n\n\n END" << std::endl;
#endif
#if 1
  auto null_resource = Register("null", std::make_unique<TraceMemoryResource>("null", null_memory_resource()));
  auto heap_resource = Register("heap", std::make_unique<TraceMemoryResource>("heap", new_delete_resource()));

  auto pool_resource = Register("pool", std::make_unique<TraceMemoryResource>("pool", std::make_unique<MemoryPoolResource>(heap_resource.memory_resource())));

  auto pinned_resource = Register("pinned", std::make_unique<TraceMemoryResource>("pinned", 
    std::make_unique<PinnedMemoryResource>(new_delete_resource())));

  host::set_default_resource(null_memory_resource());

  host::set_default_resource("heap");
  
  // constexpr std::size_t Kb = 1024;
  // constexpr std::size_t Mb = 1024*Kb;
  // std::array<std::byte,100*Kb> stack;
  // auto stack_resource = Register("stack", 
  //   std::make_unique<TraceMemoryResource>("stack", 
  //     std::make_unique<std::pmr::monotonic_buffer_resource>(std::data(stack), sizeof(stack))
  //   )
  // );

  // auto scoped_default = host::scoped_default_resource("stack");

  auto pinned_pool = Register("pinned_pool",
    std::make_unique<TraceMemoryResource>("pinned_pool",
      std::make_unique<MemoryPoolResource>(
        std::make_unique<PinnedMemoryResource>()
      )
    )
  );

  register_resource("malloc_free",
    std::make_unique<TraceMemoryResource>("malloc_free",
      std::make_unique<MemoryResourceAdaptor>(
        [](std::size_t bytes, std::size_t alignment){ return malloc(bytes); },
        [](void* ptr, std::size_t bytes, std::size_t alignment){ return free(ptr); }
      )
    )
  );

#if 1
  for( size_t j=0; j<2; ++j)
  {
    std::size_t size = 32;
    std::cout << "\n\n ITERATION " << j+1 << std::endl;
    {
      host::scoped_default_resource default_resource("pool");
      std::cout << "v1 alloc" << std::endl;
      vector<int> v1;
      v1.assign(size,1);
      //std::cout << v1.data() << " -- " << stack.data() << "   diff: " << bytes_offset(stack.data(),v1.data())  << std::endl;
      //std::cout << "    v1  = "; print_vec(std::span<int>(v1));

      //print_stacktrace(std::cout, 4, CURRENT_STACKTRACE());

      std::cout << "v2 alloc" << std::endl;
      {
        vector<int> v2;
        v2.assign(size,2);
        //std::cout << v2.data() << " -- " << stack.data() << "   diff: " << bytes_offset(stack.data(),v2.data()) << std::endl;
        //std::cout << "    v2  = "; print_vec(std::span<int>(v2));
        std::cout << "v2 dealloc" << std::endl;
      }

    // std::cout << "v3" << std::endl;
    // {
    //   vector<int> v3;
    //   v3.assign(size,3);
    //   //std::cout << v3.data() << " -- " << stack.data() << "   diff: " << bytes_offset(stack.data(),v3.data()) << std::endl;
    //   print_vec(std::span<int>(v3));
    // }

    // std::cout << "v4" << std::endl;
    // {
    //   vector<int> v4;
    //   v4.assign(size,4);
    //   //std::cout << v4.data() << " -- " << stack.data() << "   diff: " << bytes_offset(stack.data(),v4.data()) << std::endl;
    //   print_vec(std::span<int>(v4));
    // }
    // std::size_t bytes     = 48;
    // std::size_t alignment = 64;
    // auto* p = resource.allocate(bytes,alignment);
    // std::cout << p << " -- " << buffer2.data() << "   diff: " << bytes_offset(buffer2.data(),p) << std::endl;
    // resource.deallocate(p,bytes,alignment);

    if( true ) {
      //pluto::Scope pinned("pinned_pool");
      std::cout << "arr alloc" << std::endl;
      host::Array<int> arr(1000*size,256);
      arr.assign(size,5);
      //std::cout << "    arr = "; print_vec(std::span<int>(std::span<int>(arr).begin(),size));
      // std::cout << "Memory:\n";
      // pluto::TraceMemoryResource::report(std::cout);
      // std::cout << std::endl;

      // std::cout << "deallocate arr" << std::endl;
      std::cout << "arr dealloc" << std::endl;
    }

    std::cout << "v1 dealloc" << std::endl;

    // std::cout << "deallocated arr" << std::endl;
    }
    //pinned_pool.release();
    //std::cout << "released" << std::endl;
  }

  std::cout << "\n\n\n END ITERATIONS\n\n\n" << std::endl;
#endif

  // pluto::Stream stream;
  // pluto::set_default_stream(stream);

  auto managed_resource = Register("managed", std::make_unique<TraceMemoryResource>("managed", pluto::managed_resource()));
  auto device_resource = Register("device", std::make_unique<TraceMemoryResource>("device", pluto::device_resource()));

  // auto device_pool_resource = Register("device_pool", 
  //   std::make_unique<TraceMemoryResource>("device_pool",
  //     std::make_unique<std::pmr::unsynchronized_pool_resource>(&managed_resource.memory_resource()))
  // );

  auto device_pool_resource = Register("device_pool", 
      std::make_unique<TraceMemoryResource>("device_pool",
          std::make_unique<MemoryPoolResource>(device_resource)
      )
  );


  pluto::device::set_default_resource("device");
  {
    // pluto::device::scoped_default_resource default_resource("device_pool");
    // pluto::device::set_default_resource("device_pool");
    // kernel launch
    std::cout << "device arr" << std::endl;
           pluto::device::Array<int> darr(32);
  }
  pluto::wait();
  std::cout << "device arr done" << std::endl;

  // std::cout << "Memory:\n";
  // pluto::TraceMemoryResource::report(std::cout);
  // std::cout << std::endl;

  // for( int i=max_size_T; i<2*max_size_T; ++i ) {
  //   vec.push_back(i);
  // }
  // std::cout << vec.data() << " -- " << buffer.data() << std::endl;

  // std::pmr::vector<int> pmr {0,1,2};
  // const std::span std = pmr;
  //host::unregister_resource("stack");
  // host::unregister_resource("null");
  // host::unregister_resource("heap");
  
  std::cout << "\n\n\nEND\n\n\n" << std::endl;

  // std::size_t alignment = 256;
  // std::size_t aligned_down = pluto::align_down(95624137,alignment);
  // std::cout << pluto::align_down(95624137,128) << std::endl;
  // std::cout << double(aligned_down)/double(alignment) << std::endl;

  // auto compute_aligned_size = [](size_t size, size_t alignment) {
  //   size_t div           = size / alignment;
  //   size_t mod           = size % alignment;
  //   size_t _aligned_size = div * alignment;
  //   if (mod > 0) {
  //       _aligned_size += alignment;
  //   }
  //   return _aligned_size;
  // };

  // std::cout << compute_aligned_size(95624137,alignment) << std::endl; 
  // std::cout << pluto::align_up(95624137,alignment) << std::endl;

#endif
#if 1
{
  std::cout << "\n\n TEST 1\n\n" << std::endl;

  host::scoped_default_resource default_host_resource("heap");
  device::scoped_default_resource default_device_resource("device");

  host::allocator<double>   host_allocator;
  device::allocator<double> device_allocator;

  std::size_t size = 100000000;
  std::size_t bytes = size * sizeof(double);
  double* h_ptr = host_allocator.allocate(size);
  double* d_ptr = device_allocator.allocate(size);

  std::cout << "h_ptr = " << h_ptr << std::endl;
  std::cout << "d_ptr = " << d_ptr << std::endl;

  std::cout << "is_host(h_ptr)    : " << is_host(h_ptr) << std::endl;
  std::cout << "is_host(d_ptr)    : " << is_host(d_ptr) << std::endl;
  std::cout << "is_device(h_ptr)  : " << is_device(h_ptr) << std::endl;
  std::cout << "is_device(d_ptr)  : " << is_device(d_ptr) << std::endl;
  std::cout << "is_managed(h_ptr) : " << is_managed(h_ptr) << std::endl;
  std::cout << "is_managed(d_ptr) : " << is_managed(d_ptr) << std::endl;

  std::fill(h_ptr, h_ptr+size, 1.);
  // auto& stream = get_default_stream();
  Stream stream1, stream2, stream3;
  memcpy_host_to_device(d_ptr, h_ptr, bytes, stream1);

  set_on_device(stream2, d_ptr, size, 2.);

  memcpy_device_to_host(h_ptr, d_ptr, bytes, stream3);
  stream3.wait();

  std::cout << "*********** h_ptr[size-1] = " << h_ptr[size-1] << std::endl;

  host_allocator.deallocate(h_ptr, size);
  device_allocator.deallocate(d_ptr, size);
}

if (0) {
  std::cout << "\n\n TEST 2: managed\n\n" << std::endl;
  host::scoped_default_resource default_resource("managed");
  host::allocator<double>   host_allocator;
  // device::allocator<double> device_allocator;

  std::size_t size = 1024;
  // std::size_t bytes = size * sizeof(double);
  double* h_ptr = host_allocator.allocate(size);
  // double* d_ptr = device_allocator.allocate(size);

  std::cout << "is_host(h_ptr)    : " << is_host(h_ptr) << std::endl;
  std::cout << "is_device(h_ptr)  : " << is_device(h_ptr) << std::endl;
  std::cout << "is_managed(h_ptr) : " << is_managed(h_ptr) << std::endl;

  std::fill(h_ptr, h_ptr+size, 1.);

  // memcpy_host_to_device(d_ptr, h_ptr, bytes);

  set_on_device(h_ptr, size, 2.);
  pluto::wait();

  // memcpy_device_to_host(h_ptr, d_ptr, bytes);

  std::cout << "h_ptr[size-1] = " << h_ptr[size-1] << std::endl;

  host_allocator.deallocate(h_ptr, size);
  // device_allocator.deallocate(d_ptr, size);
}

if (0)
{
  std::cout << "\n\n TEST 3: pinned+mapped\n\n" << std::endl;
  host::scoped_default_resource default_resource("pinned_pool");

  host::allocator<double>   host_allocator;
  std::size_t size = 1024;
  std::size_t bytes = size * sizeof(double);
  double* h_ptr = host_allocator.allocate(size);
  double* d_ptr = get_registered_device_pointer(h_ptr);

  std::cout << "h_ptr = " << h_ptr << std::endl;
  std::cout << "d_ptr = " << d_ptr << std::endl;

  std::cout << "is_host(h_ptr)    : " << is_host(h_ptr) << std::endl;
  std::cout << "is_host(d_ptr)    : " << is_host(d_ptr) << std::endl;
  std::cout << "is_device(h_ptr)  : " << is_device(h_ptr) << std::endl;
  std::cout << "is_device(d_ptr)  : " << is_device(d_ptr) << std::endl;
  std::cout << "is_managed(h_ptr) : " << is_managed(h_ptr) << std::endl;
  std::cout << "is_managed(d_ptr) : " << is_managed(d_ptr) << std::endl;

  std::fill(h_ptr, h_ptr+size, 1.);

  memcpy_host_to_device(d_ptr, h_ptr, bytes);
  pluto::wait();

  // launch_kernel(size, [=] HIC_DEVICE (auto i) {
  //   d_ptr[i] = i;
  // });
  //
  pluto::wait();

  memcpy_device_to_host(h_ptr, d_ptr, bytes);

  std::cout << "h_ptr[size-1] = " << h_ptr[size-1] << std::endl;

  host_allocator.deallocate(h_ptr, size);
  // device_allocator.deallocate(d_ptr, size);
}

if(0)
{
  std::cout << "\n\n TEST 4: pinned nonmapped\n\n" << std::endl;
  host::scoped_default_resource default_host_resource("pinned_pool");
  device::scoped_default_resource default_device_resource("device");
  host::allocator<double> host_allocator;
  device::allocator<double> device_allocator;

  std::size_t size = 1025;
  std::size_t bytes = size * sizeof(double);
  double* h_ptr = host_allocator.allocate(size);
  double* d_ptr = device_allocator.allocate(size);

  std::cout << "h_ptr = " << h_ptr << std::endl;
  std::cout << "d_ptr = " << d_ptr << std::endl;

  std::cout << "is_host(h_ptr)    : " << is_host(h_ptr) << std::endl;
  std::cout << "is_host(d_ptr)    : " << is_host(d_ptr) << std::endl;
  std::cout << "is_device(h_ptr)  : " << is_device(h_ptr) << std::endl;
  std::cout << "is_device(d_ptr)  : " << is_device(d_ptr) << std::endl;
  std::cout << "is_managed(h_ptr) : " << is_managed(h_ptr) << std::endl;
  std::cout << "is_managed(d_ptr) : " << is_managed(d_ptr) << std::endl;

  std::fill(h_ptr, h_ptr+size, 1.);

  memcpy_host_to_device(d_ptr, h_ptr, bytes);

  // launch_kernel(size, [=] HIC_DEVICE (auto i) {
  //   d_ptr[i] = i;
  // });

  pluto::wait();

  memcpy_device_to_host(h_ptr, d_ptr, bytes);

  std::cout << "h_ptr[size-1] = " << h_ptr[size-1] << std::endl;

  host_allocator.deallocate(h_ptr, size);
  device_allocator.deallocate(d_ptr, size);
}

#endif

for( int i=0; i<2; ++i)
{
  std::cout << "\n\n TEST 5 device memory pool ITERATION " << i+1 << "\n\n" << std::endl;

  host::scoped_default_resource default_host_resource("heap");
  device::scoped_default_resource default_device_resource("device_pool");

  host::allocator<double>   host_allocator;
  device::allocator<double> device_allocator;

  std::size_t size = 1024;
  double* h_ptr = host_allocator.allocate(size);
  double* d_ptr = device_allocator.allocate(size);

  std::cout << "h_ptr = " << h_ptr << std::endl;
  std::cout << "d_ptr = " << d_ptr << std::endl;

  std::cout << "is_host(h_ptr)    : " << is_host(h_ptr) << std::endl;
  std::cout << "is_host(d_ptr)    : " << is_host(d_ptr) << std::endl;
  std::cout << "is_device(h_ptr)  : " << is_device(h_ptr) << std::endl;
  std::cout << "is_device(d_ptr)  : " << is_device(d_ptr) << std::endl;
  std::cout << "is_managed(h_ptr) : " << is_managed(h_ptr) << std::endl;
  std::cout << "is_managed(d_ptr) : " << is_managed(d_ptr) << std::endl;

  std::fill(h_ptr, h_ptr+size, 1.);

  copy_host_to_device(d_ptr, h_ptr, size);

  set_on_device(d_ptr, size, 2.);

  copy_device_to_host(h_ptr, d_ptr, size);

  std::cout << "*********** h_ptr[size-1] = " << h_ptr[size-1] << std::endl;

  host_allocator.deallocate(h_ptr, size);
  device_allocator.deallocate(d_ptr, size);
  std::cout << "\n\n END TEST 5 device memory pool\n\n" << std::endl;
}

#define PLUTO_KERNEL [=] HIC_DEVICE

// Custom type
{
  std::cout << "\n\n TEST 6 custom type\n\n" << std::endl;
  device::scoped_default_resource default_device_resource("device_pool");
  
  bool verbosity = TraceOptions::instance().enabled;

  {
    std::cout << "Using device::allocator<CustomType>" << std::endl;
    auto alloc = device::allocator<CustomType>{};
    CustomType* x = alloc.allocate(1);
    alloc.construct(x,true);
    launch_kernel( PLUTO_KERNEL { x->compute(); });
    alloc.destroy(x);
    alloc.deallocate(x,1);
  }
  {
    std::cout << "\n\ndevice::make_unique<CustomType>" << std::endl;
    auto unique_x = device::make_unique<CustomType>(verbosity);
    auto x = unique_x.get(); 
    launch_kernel(PLUTO_KERNEL { x->compute(); });
    CustomType ht;
    copy_device_to_host(&ht, x);
   }
   { 
     std::cout << "\n\nhost::make_copy, device::make_copy" << std::endl;
     auto dt  = device::make_unique<CustomType>(verbosity);
     auto ht  = host::make_copy(dt);
     ht->a = 20;
     auto dt2 = device::make_copy(ht);
  }
   {
    std::cout << "\n\ncopy_host_to_device, copy_device_to_host" << std::endl;
	   auto ht = host::make_unique<double>();
	   *ht = 9.4;
	   auto dt = device::make_copy(ht);
	   auto ht2 = host::make_copy(dt);
     assert(*ht2 == 9.4);
     *ht2 = 0.3;
     copy_host_to_device(dt,ht2);
     copy_device_to_host(ht,dt);
     assert(*ht == 0.3);
   }
}
  std::cout << "\n\n END TEST 6 custom type\n\n" << std::endl;

{
  std::cout << "\n\n TEST 7 device vector\n\n" << std::endl;
  device::vector<CustomType> dvec(device_resource.memory_resource());
  dvec.resize(2);
  auto ddata = dvec.data();
  launch_kernel(dvec.size(), PLUTO_KERNEL (int i) {
    ddata[i].a = i;
  });

  host::vector<CustomType> hvec(heap_resource.memory_resource());
  hvec.resize(4);
  copy_device_to_host(hvec.data(), dvec.data(), dvec.size());
  assert( hvec[0].a == 0 );
  assert( hvec[1].a == 1 );
  assert( hvec[2].a == 999 );
  assert( hvec[3].a == 999 );
}
 std::cout << "\n\n END TEST 7 device vector\n\n" << std::endl;

return 0;
}
