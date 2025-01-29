

#include <vector>
//#include <span>
#include <iostream>
#include <memory>
#include <new>
#include <cstdlib>
#include <map>

#include "pluto/pluto.h"
#include "benchmark.h"

static std::size_t to_bytes(const std::string& str) {
    auto unit_to_bytes = [](char unit) {
       static const std::map<char, std::size_t> map{
            {'G',1024*1024*1024},
            {'M',1024*1024},
            {'K',1024},
            {'B',1}
        };
        return map.at(static_cast<char>(std::toupper(unit)));
    };
    for (char unit: {'G','g','M','m','K','k','B','b'}) {
        if (auto pos = str.find(unit); pos != std::string::npos) {
            return unit_to_bytes(unit) * std::stoull(str.substr(0,pos));
        }
    }
    return std::stoull(str);
}

void setup_resources(std::size_t bytes) {
  using namespace pluto;
  auto null_resource   = register_resource("null",
    std::make_unique<TraceMemoryResource>("null", null_memory_resource())
  );
  auto heap_resource   = register_resource("heap",
    std::make_unique<TraceMemoryResource>("heap", new_delete_resource())
  );
  auto pinned_resource = register_resource("pinned",
    std::make_unique<TraceMemoryResource>("pinned", 
      std::make_unique<PinnedMemoryResource>(new_delete_resource())
    )
  );
  auto pinned_pool_resource = register_resource("pinned_pool",
    std::make_unique<TraceMemoryResource>("pinned_pool",
      std::make_unique<MemoryPoolResource>(
        std::make_unique<TraceMemoryResource>("pinned", 
          std::make_unique<PinnedMemoryResource>(new_delete_resource())
        )
      )
    )
  );
  pinned_pool_resource->deallocate(pinned_pool_resource->allocate(4*bytes), 4*bytes);

  auto managed_resource = register_resource("managed",
    std::make_unique<TraceMemoryResource>("managed", pluto::managed_resource())
  );
  auto device_resource = register_resource("device",
    std::make_unique<TraceMemoryResource>("device", pluto::device_resource())
  );

  auto device_pool_resource = register_resource("device_pool", 
      std::make_unique<TraceMemoryResource>("device_pool",
        std::make_unique<MemoryPoolResource>(device_resource)
      )
  );
  device_pool_resource->deallocate(device_pool_resource->allocate(4*bytes), 4*bytes);
}
bool TRACE() {
    char* val;                                                                        
    val = getenv( "TRACE" );                                                       
    int retval = 0;
    if (val != NULL) {                                                                 
        retval = std::atoi(val);
    }
    return retval;                                                                        
}

int main(int argc, char* argv[]) {
  using namespace pluto;
  using namespace pluto;
  using value_type = double;

  TraceOptions::instance().enabled = TRACE();
    pluto::set_trace(TRACE());

  std::size_t bytes = 1024*1024; // 1 Mb
  if( argc > 1 ) {
    bytes = to_bytes(argv[1]);
  }
  int n = static_cast<int>(bytes/sizeof(value_type));
  int num_repeats = 25;
  int num_warmups = 5;
  if( argc > 2 ) {
    num_repeats = std::atoi(argv[2]);
  }
  if( argc > 3 ) {
    num_warmups = std::atoi(argv[3]);
  }
  setup_resources(bytes);

{
    value_type const v_input_1{1.0};
    value_type const v_input_2{1.0};
    value_type const v_output{0.0};
    value_type const v_output_reference{v_input_1 + v_input_2};

    double d2h_volume  = double(1 * n * sizeof(value_type)) / 1024./1024.; // Mb
    double h2d_volume  = double(2 * n * sizeof(value_type)) / 1024./1024.; // Mb
    double copy_volume = double(3 * n * sizeof(value_type)) / 1024./1024.; // Mb
    std::cout << "copy_volume = " << copy_volume << " Mb" << std::endl;

    value_type* h_input_1;
    value_type* h_input_2;
    value_type* h_output;
    value_type* d_input_1;
    value_type* d_input_2;
    value_type* d_output;

    auto allocate_host = [&] {
      auto start = std::chrono::steady_clock::now();
      host::allocator<value_type>   allocator;
      h_input_1 = allocator.allocate(n);
      h_input_2 = allocator.allocate(n);
      h_output  = allocator.allocate(n);
      auto end = std::chrono::steady_clock::now();
      std::cout << "\tallocate_host    = " << std::chrono::duration<float>(end-start).count() << std::endl;
    };
    auto allocate_device = [&] {
      auto start = std::chrono::steady_clock::now();
      device::allocator<value_type> allocator;
      d_input_1 = allocator.allocate(n);
      d_input_2 = allocator.allocate(n);
      d_output  = allocator.allocate(n);
      auto end = std::chrono::steady_clock::now();
      std::cout << "\tallocate_device  = " << std::chrono::duration<float>(end-start).count() << std::endl;
    };
    auto use_registered_device_pointers = [&] {
      d_input_1 = pluto::get_registered_device_pointer(h_input_1);
      d_input_2 = pluto::get_registered_device_pointer(h_input_2);
      d_output  = pluto::get_registered_device_pointer(h_output);
    };
    auto deallocate_host = [&] {
      auto start = std::chrono::steady_clock::now();
      host::allocator<value_type>   allocator;
      allocator.deallocate(h_input_1, n);
      allocator.deallocate(h_input_2, n);
      allocator.deallocate(h_output, n);
      auto end = std::chrono::steady_clock::now();
      std::cout << "\tdeallocate_host  = " << std::chrono::duration<float>(end-start).count() << std::endl;
    };
    auto deallocate_device = [&] {
      auto start = std::chrono::steady_clock::now();
      device::allocator<value_type> allocator;
      allocator.deallocate(d_input_1, n);
      allocator.deallocate(d_input_2, n);
      allocator.deallocate(d_output, n);
      auto end = std::chrono::steady_clock::now();
      std::cout << "\tdeallocate_device = " << std::chrono::duration<float>(end-start).count() << std::endl;
    };
    auto initialize_host = [&]() {
      initialize_host_memory(h_input_1, n, v_input_1);
      initialize_host_memory(h_input_2, n, v_input_2);
      initialize_host_memory(h_output,  n, v_output);
    };

    float ref_time = 0.;
    float ref_host_time = 0.;
    float ref_offload_time = 0.;
    float ref_device_time = 0.;
    float ref_copy_time = 0.;
    float ref_d2h_time = 0.;
    float ref_h2d_time = 0.;
    float time = 0.;
    float offload_time = 0.;
    float host_time = 0;
    float device_time = 0;
    float copy_time = 0;
    float d2h_time = 0;
    float h2d_time = 0;
    int count = 0;
    int count_time = 0;
    auto timers_init = [&]{
      host_time = 0;
      device_time = 0;
      copy_time = 0;
      d2h_time = 0;
      h2d_time = 0;
      count = 0;
      count_time = 0;
    };
    auto timers_final = [&]{
      host_time /= float(count_time);
      device_time /= float(count_time);
      d2h_time /= float(count_time);
      h2d_time /= float(count_time);
      copy_time = d2h_time + h2d_time;
      if (!ref_time) {
        ref_time = time;
      }
      if (!ref_host_time) {
        ref_host_time = host_time;
      }
      if (!ref_device_time) {
        ref_device_time = device_time;
      }
      if (!ref_copy_time) {
        ref_copy_time = copy_time;
      }
      if (!ref_d2h_time) {
        ref_d2h_time = d2h_time;
      }
      if (!ref_h2d_time) {
        ref_h2d_time = h2d_time;
      }


      offload_time = time - host_time;
      // if(host_time/ref_host_time > 1.1) {
        offload_time = time - ref_host_time;
      // }
      if (!ref_offload_time) {
        ref_offload_time = offload_time;
      }

      auto unaccounted_time = time - host_time - device_time - copy_time;
      std::cout << "\ttotal_time       = " << time << " ( " << time/ref_time << " ref ) " << std::endl;
      std::cout << "\toffload_time     = " << offload_time << " ( " << offload_time/ref_offload_time << " ref ) " << std::endl;
      std::cout << "\thost_time        = " << host_time << " ( " << host_time/time << " total_time , " << host_time/ref_host_time << " ref )" <<  std::endl;
      std::cout << "\tdevice_time      = " << device_time << " ( " << device_time/offload_time << " offload_time , " << device_time/ref_device_time << " ref ) " <<  std::endl;
      if (copy_time == 0. || (ref_device_time != 0 && device_time/ref_device_time > 1.1) ) {
        auto impl_d2h_time = d2h_time + (host_time - ref_host_time);
        auto impl_h2d_time = h2d_time + (device_time - ref_device_time);
        auto impl_copy     = impl_d2h_time + impl_h2d_time;

        std::cout << "\timpl_d2h_time    = " << impl_d2h_time << " ( " << impl_d2h_time/offload_time << " offload_time , " << impl_d2h_time/ref_d2h_time << " ref, "  << d2h_volume/impl_d2h_time/1024. << " Gb/s ) " <<  std::endl;
        std::cout << "\timpl_h2d_time    = " << impl_h2d_time << " ( " << impl_h2d_time/offload_time << " offload_time , " << impl_h2d_time/ref_h2d_time << " ref, "  << h2d_volume/impl_h2d_time/1024. << " Gb/s ) " <<  std::endl;
        std::cout << "\timpl_copy_time   = " << impl_copy     << " ( " << impl_copy/offload_time     << " offload_time , " << impl_copy/ref_copy_time    << " ref, " << copy_volume/impl_copy/1024.    << " Gb/s ) " <<  std::endl;
      }
      else {
        std::cout << "\td2h_time         = " << d2h_time   << " ( " << d2h_time/offload_time  << " offload_time , " << d2h_time/ref_d2h_time   << " ref, "  << d2h_volume/d2h_time/1024.   << " Gb/s ) " <<  std::endl;
        std::cout << "\th2d_time         = " << h2d_time   << " ( " << h2d_time/offload_time  << " offload_time , " << h2d_time/ref_h2d_time   << " ref, "  << h2d_volume/h2d_time/1024.   << " Gb/s ) " <<  std::endl;
        std::cout << "\tcopy_time        = " << copy_time  << " ( " << copy_time/offload_time << " offload_time , " << copy_time/ref_copy_time << " ref, " << copy_volume/copy_time/1024. << " Gb/s ) " <<  std::endl;
      }
      std::cout << "\tunaccounted_time = " << unaccounted_time << " ( " << unaccounted_time/offload_time << " offload_time )" << std::endl;

    };

    bool in_verify = false;
    auto touch_host = [&](value_type* host_buffer) -> value_type {
      if (not in_verify) {
        return touch_host_memory(host_buffer, n);
      }
      return 0.;
    };

    auto verify = [&](std::function<void()> f) {
      in_verify = true;
      initialize_host();
      f();
      if (verify_host_memory(h_output, n, v_output_reference) == false) {
         std::cout << "\tVerification FAILED\n" << std::endl;
      }
      initialize_host();
      in_verify = false;
    };

    auto measure_performance = [&](std::function<void()> f) {
      timers_init();
      time = run_timed(f,num_repeats,num_warmups);
      timers_final();
    };

{
    std::cout << "\n BENCHMARK unpinned (reference)\n" << std::endl;

    host::scoped_default_resource   default_host_resource("heap");
    device::scoped_default_resource default_device_resource("device");
    allocate_host();
    allocate_device();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};
      float d2h{0};
      float h2d{0};
      htime += touch_host(h_input_1);
      htime += touch_host(h_input_2);
      h2d   += run_timed([&]{copy_host_to_device(d_input_1, h_input_1, n);});
      h2d   += run_timed([&]{copy_host_to_device(d_input_2, h_input_2, n);});
      dtime += launch_benchmark(d_output, d_input_1, d_input_2, n);
      d2h   += run_timed([&]{copy_device_to_host(h_output, d_output, n);});
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time   += htime;
          device_time += dtime;
          d2h_time    += d2h;
          h2d_time    += h2d;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    deallocate_device();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}


{
    std::cout << "\n BENCHMARK pinned\n" << std::endl;

    host::scoped_default_resource   default_host_resource("pinned");
    device::scoped_default_resource default_device_resource("device");
    allocate_host();
    allocate_device();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};
      float d2h{0};
      float h2d{0};
      htime += touch_host(h_input_1);
      htime += touch_host(h_input_2);
      h2d   += run_timed([&]{copy_host_to_device(d_input_1, h_input_1, n);});
      h2d   += run_timed([&]{copy_host_to_device(d_input_2, h_input_2, n);});
      dtime += launch_benchmark(d_output, d_input_1, d_input_2, n);
      d2h   += run_timed([&]{copy_device_to_host(h_output, d_output, n);});
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time   += htime;
          device_time += dtime;
          d2h_time    += d2h;
          h2d_time    += h2d;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    deallocate_device();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}

{
    std::cout << "\n BENCHMARK pinned_pool + device_pool \n" << std::endl;
    host::scoped_default_resource   default_host_resource("pinned_pool");
    device::scoped_default_resource default_device_resource("device_pool");
    allocate_host();
    allocate_device();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};
      float d2h{0};
      float h2d{0};
      htime += touch_host(h_input_1);
      htime += touch_host(h_input_2);
      h2d   += run_timed([&]{copy_host_to_device(d_input_1, h_input_1, n);});
      h2d   += run_timed([&]{copy_host_to_device(d_input_2, h_input_2, n);});
      dtime += launch_benchmark(d_output, d_input_1, d_input_2, n);
      d2h   += run_timed([&]{copy_device_to_host(h_output, d_output, n);});
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time   += htime;
          device_time += dtime;
          d2h_time    += d2h;
          h2d_time    += h2d;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    deallocate_device();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}


{
    std::cout << "\n BENCHMARK managed\n" << std::endl;

    host::scoped_default_resource   default_host_resource("managed");
    allocate_host();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};

      htime += touch_host(h_input_1);
      htime += touch_host(h_input_2);
      dtime += launch_benchmark(h_output, h_input_1, h_input_2, n);
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time += htime;
          device_time += dtime;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}

{
    std::cout << "\n BENCHMARK managed with prefetch\n" << std::endl;

    host::scoped_default_resource   default_host_resource("managed");
    allocate_host();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};
      float h2d{0};
      float d2h{0};
      htime += touch_host(h_input_1);
      h2d   += run_timed([&]{prefetch_host_to_device(h_input_1, n*sizeof(value_type));});
      htime += touch_host(h_input_2);
      h2d   += run_timed([&]{prefetch_host_to_device(h_input_2, n*sizeof(value_type));});
      h2d   += run_timed([&]{pluto::wait();});
      dtime += launch_benchmark(h_output, h_input_1, h_input_2, n);
      d2h   += run_timed([&]{
        prefetch_device_to_host(h_output, n*sizeof(value_type));
        pluto::wait();
      });
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time += htime;
          device_time += dtime;
          d2h_time    += d2h;
          h2d_time    += h2d;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}


{
    std::cout << "\n BENCHMARK zero-copy\n" << std::endl;

    host::scoped_default_resource   default_host_resource("pinned");
    allocate_host();
    use_registered_device_pointers();
    initialize_host();

    auto function = [&]() {
      float htime{0};
      float dtime{0};

      htime += touch_host(h_input_1);
      htime += touch_host(h_input_2);
      dtime += launch_benchmark(d_output, d_input_1, d_input_2, n);
      htime += touch_host(h_output);
      if( count >= num_warmups ) {
          ++count_time;
          host_time += htime;
          device_time += dtime;
      }
      ++count;
    };

    verify(function);
    measure_performance(function);

    deallocate_host();
    std::cout << "\n END BENCHMARK \n" << std::endl;
}
}

  unregister_resources();

return 0;
}
