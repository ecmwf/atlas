
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "kernel.h"

#include "pluto/memory_resource.h"
#include "hic/hic.h"



// template<typename Base, typename Derived, typename... Args>
// void construct_on_device(Base** ptr, Args... args) {
//     launch_kernel(1, [=] HIC_DEVICE(int i) {
//         // printf("arg1 = %d\n", args);
//         *ptr = new Derived(args...);
//     });
// }

template <typename T>
class DeviceObject;

template <typename T>
class DevicePtr {
public:
    HIC_HOST_DEVICE
    DevicePtr() {
        // device_ptr = nullptr;
// #ifndef __CUDA_ARCH__
//         std::cout << "DevicePtr() " <<  this << std::endl;
// #endif
    }

//     HICs_HOST_DEVICE
//     DevicePtr(const DevicePtr& other) {
// // #ifndef __CUDA_ARCH__
// //         std::cout << "DevicePtr(other) " <<  this << std::endl;
// // #endif
// //         device_ptr = other.device_ptr;
//     }

	template< typename T1, typename = std::enable_if_t<std::is_convertible_v<T1*, T*>> >
	HIC_HOST_DEVICE DevicePtr(const DevicePtr<T1> &dp)
		: device_ptr((T1*)dp)
	{}

	template< typename T1, typename = std::enable_if_t<std::is_convertible_v<T1*, T*>> >
	HIC_HOST DevicePtr(const DeviceObject<T1> &dp)
		: device_ptr((T1*)dp)
	{}

//      HIC_HOST_DEVICE
//      ~DevicePtr() {
// // #ifndef __CUDA_ARCH__
// // #warning host
// //         std::cout << "~DevicePtr() " <<  this << std::endl;
// //         //clear();
// //         //cudaFree(device_ptr);
// // #else
// // #warning device
// // #endif
//     }


    // template<typename Derived, typename... Args>
    // HIC_HOST
    // void construct(Args... args) {
    //     //cudaMalloc(&device_ptr, sizeof(T*));
    //     launch_construct<Derived>(device_ptr, std::forward<Args>(args)...);
    // }

    // HIC_HOST
    // void clear() {
    //     std::cout << "  clear " <<  this << std::endl;
    //     if( !device_ptr ) {
    //         std::cout << "device_ptr already null" << std::endl;
    //     }
    //     else {
    //         launch_destruct(device_ptr);
    //         // cudaFree(device_ptr);
    //         device_ptr = nullptr;
    //     }
    // }


    HIC_HOST_DEVICE
    T* get() const { return device_ptr; }

    HIC_HOST_DEVICE
    T* operator ->() const { return device_ptr; }

    HIC_HOST_DEVICE
    T& operator *() const { return *device_ptr; }

	HIC_HOST_DEVICE operator T*() const { return  device_ptr; }


    // HIC_HOST_DEVICE
    // T** double_ptr() const { return device_ptr; }


    HIC_HOST static DevicePtr FromRawDevicePtr(T *p) {
        return DevicePtr{ p };
    }
private:
    HIC_HOST_DEVICE __inline__ explicit DevicePtr(T *p) : device_ptr(p) {}
    T* device_ptr = nullptr;
};


template<typename T>
class DeviceMemory {
    T *_p = nullptr;
    std::size_t _bytes;
    DeviceMemory(std::size_t bytes) : _bytes(bytes) { cudaMalloc(&_p, _bytes); }
public:
    static DeviceMemory AllocateElements(std::size_t n) {return {n*sizeof(T)}; }
    static DeviceMemory AllocateBytes(std::size_t bytes) {return {bytes}; }
    ~DeviceMemory() { if (_p) {cudaFree(_p);} }

    HIC_HOST
    operator DevicePtr<T>() const {
        return DevicePtr<T>::FromRawDevicePtr(_p);
    }
    operator bool() const {
        return (_p != nullptr);
    }
};

template<typename T, typename... Args> 
HIC_HOST_DEVICE void allocate_object(DevicePtr<T> p, Args... args) {
    new (p.get()) T(args...);
}

template<typename T, typename... Args> 
HIC_HOST_DEVICE void new_on_device(T* p, Args... args) {
    new (p) T(args...);
}

template<typename T>
HIC_HOST_DEVICE void delete_object(DevicePtr<T> p) {
    p->~T();
}

template<typename T>
HIC_HOST_DEVICE void delete_on_device(T* p) {
    p->~T();
}


template<typename T>
class DeviceObject {
public:
    T* ptr_ {nullptr};
public:
    template<typename... Args>
    DeviceObject(Args... args) { 
        std::cout << "allocate object on device" << std::endl;
#if HIC_COMPILER
        HIC_CALL( hicMalloc(&ptr_, sizeof(T)) );
        new_on_device<T><<<1, 1>>>(ptr_, args...);
        HIC_CALL( hicDeviceSynchronize() );
#else
        ptr_ = (T*)malloc(sizeof(T)); new_on_device<T>(ptr_, args...);
#endif

    }

     DeviceObject& operator=( const DeviceObject& ) = delete; // non copyable

    ~DeviceObject() {
        if (ptr_) {
            std::cout << "deallocate object on device" << std::endl;
#if HIC_COMPILER
            delete_on_device<T><<<1, 1>>>(ptr_);
            HIC_CALL( hicDeviceSynchronize() );
            HIC_CALL( hicFree(ptr_) );
#else
            delete_on_device<T>(ptr_);
            free(ptr_);
#endif
        }
    }

    HIC_HOST
    operator DevicePtr<T>() const { 
        return DevicePtr<T>::FromRawDevicePtr(ptr_);
    }

    HIC_HOST
    operator T*() const { 
        return ptr_;
    }

};

HIC_DEVICE pluto::memory_resource* gpu_default_;

// __global__ void set_default_global(pluto::memory_resource* ptr) {
//     gpu_default_= DevicePtr<pluto::memory_resource>::FromRawDevicePtr(ptr);
// }

HIC_DEVICE void set_device_ptr(pluto::memory_resource*& destination, pluto::memory_resource* source) {
    destination = source;
}

#if HIC_COMPILER
__global__ void launch_set_default(pluto::memory_resource* ptr) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0) {
        set_device_ptr(gpu_default_, ptr);
    }
}
#endif

void set_default_on_device(pluto::memory_resource* ptr) {
#if HIC_COMPILER
    launch_set_default<<<1,1>>>(ptr);
    HIC_CHECK_KERNEL_LAUNCH();
#else
    gpu_default_ = ptr;
#endif
}

void set_on_device(double* x, std::size_t size, double value) {



    // auto& resource = pluto::device::new_delete_resource();
    // pluto::device::allocator_on_device<double> allocator(resource);


    // pluto::device::new_delete_resource_t resource;

    // using memory_resource = pluto::memory_resource;
    // using new_delete_resource_t = pluto::device::new_delete_resource_t;
    // resource_t* host_mr;
    // host_mr = new resource_t();

    // Option 1
    // memory_resource** device_mr = nullptr;
    // cudaMalloc(&device_mr, sizeof(memory_resource*));
    // launch_kernel(1, [=] HIC_DEVICE(auto i) {
    //     *device_mr = new new_delete_resource_t();
    // });


    // Option 2
    // memory_resource** device_mr = nullptr;
    // device_mr = new (memory_resource*);
    // cudaHostRegister(&device_mr, sizeof(memory_resource*), cudaHostRegisterMapped);
    // *device_mr = new new_delete_resource_t();
    // cudaHostRegister(device_mr, sizeof(new_delete_resource_t), cudaHostRegisterMapped);

    // Option 3
    // pluto::device::cudaMallocMemoryResource from_host;
    // GPUClonable<pluto::device::cudaMallocMemoryResource> gpu_new_delete( &from_host );
    // gpu_new_delete.updateDevice();


    // Option 4
    // DevicePtr<memory_resource> device_mr;
    // device_mr.construct<new_delete_resource_t>();

    // DevicePtr<memory_resource> traced_device_mr;
    // traced_device_mr.construct<pluto::device::traced_memory_resource>(device_mr.double_ptr());

    // DeviceObject<new_delete_resource_t> obj_resource;
    // DevicePtr<memory_resource> ptr(obj_resource);
    // set_default_on_device(ptr.get());

    // gpu_default_ = ptr.get();
    // DevicePtr<memory_resource> gpu_default = DevicePtr<memory_resource>::FromRawDevicePtr(gpu_default_);

    // gpu_default_ = DevicePtr<memory_resource>(obj_resource);
    // DevicePtr<memory_resource> abstract(ptr);
    // DevicePtr<memory_resource> abstract(obj_resource);

    // std::cout << "resource initialized " << std::endl;

    launch_kernel(size, [=] HIC_DEVICE(auto i) {
        x[i] = value;
    });

    // device_mr.clear();
}

void set_on_device(const pluto::stream& stream, double* x, std::size_t size, double value) {
    std::cout << "set_on_device async" << std::endl;
    launch_kernel(stream, size, [=] HIC_DEVICE(auto i) {
        x[i] = value;
    });

    // device_mr.clear();
}
