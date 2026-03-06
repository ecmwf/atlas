/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MPIMemoryResource.h"

#include <iostream>
#include <string_view>
#include <cstdlib>

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/memory.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

#if PLUTO_HAVE_MPI
#include <mpi.h>
#endif

namespace pluto {

namespace {
    bool mpi_execution_detected() {
        static bool detected = [] {
            // Check common MPI environment variables to detect if we are running in an MPI environment.
            std::vector<std::string> mpi_detection_env_vars{
                "OMPI_COMM_WORLD_SIZE",   // OpenMPI
                "ALPS_APP_PE",            // Cray aprun
                "PMI_SIZE",               // Intel MPI
                "SLURM_STEP_NUM_TASKS"    // slurm srun
            };
            for (const auto& env : mpi_detection_env_vars) {
                if (std::getenv(env.c_str())) {
                    std::cerr << "mpi_execution_detected = true" << std::endl;
                    return true;
                }
            }
            return false;
        }();
        return detected;
    }
    bool use_mpi() {
        static bool detected = [] {
            #if PLUTO_HAVE_MPI
                int is_mpi_initialized = 0;
                MPI_Initialized(&is_mpi_initialized);
                if (is_mpi_initialized) {
                    return true;
                }
            #endif
            if (mpi_execution_detected()) {
                pluto::trace::out <<
                    "PLUTO_WARNING: MPI environment detected but MPI does not seem to be initialized.\n"
                    "               Please initialize MPI before using pluto::mpi_resource.\n"
                    "               Continuing without MPI_Alloc_mem / MPI_Free_mem support."
                    << std::endl;
            }
            return false;
        }();
        return detected;
    }
    void* mpi_allocate([[maybe_unused]] std::size_t bytes) {
        void* ptr = nullptr;
        #if PLUTO_HAVE_MPI
            int err = MPI_Alloc_mem(bytes, MPI_INFO_NULL, &ptr); 
            if (err != MPI_SUCCESS) {
                throw std::runtime_error("MPI_Alloc_mem failed");
            }
        #endif
        return ptr;
    }
    void mpi_deallocate([[maybe_unused]] void* ptr) {
        #if PLUTO_HAVE_MPI
            int err = MPI_Free_mem(ptr);
            if (err != MPI_SUCCESS) {
                throw std::runtime_error("MPI_Free_mem failed");
            }
        #endif
    }
}

void mpi_init() {
    if (mpi_execution_detected()) {
        #if PLUTO_HAVE_MPI
            int mpi_err = MPI_Init(nullptr, nullptr);
            if (mpi_err != MPI_SUCCESS) {
                throw std::runtime_error("MPI_Init failed");
            }
        #endif
    }
}

void mpi_finalize() {
    if (mpi_execution_detected()) {
        #if PLUTO_HAVE_MPI
            int mpi_err = MPI_Finalize();
            if (mpi_err != MPI_SUCCESS) {
                throw std::runtime_error("MPI_Finalize failed");
            }
        #endif
    }
}

// --------------------------------------------------------------------------------------------------------

class MPIMemoryResource : public memory_resource {
public:
    MPIMemoryResource() = default;

    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;

    bool do_is_equal(const memory_resource& other) const noexcept override;
};

// --------------------------------------------------------------------------------------------------------

void* MPIMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    void* ptr;
    if (use_mpi()) {
        ptr = mpi_allocate(bytes);
    }
    else {
        ptr = new_delete_resource()->allocate(bytes, alignment);
    }

    memory::mpi.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate(get_label(), ptr, bytes, alignment, "pluto::mpi_resource", &memory::mpi);
    }
    return ptr;
}

void MPIMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    memory::mpi.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate(get_label(), ptr, bytes, alignment, "pluto::mpi_resource", &memory::mpi);
    }

    if (use_mpi()) {
        mpi_deallocate(ptr);
    }
    else {
        new_delete_resource()->deallocate(ptr, bytes, alignment);
    }
}

bool MPIMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

// --------------------------------------------------------------------------------------------------------

namespace {
template<typename T>
struct constant_init {
    union {
        T obj;
    };
    constexpr constant_init() : obj() { }

    template<typename ...Args>
    explicit constexpr constant_init(Args... args) : obj(args...) { }

    ~constant_init() { /* do nothing, union member is not destroyed */ }
};
}

constant_init<MPIMemoryResource>  mpi_res{};

memory_resource* mpi_resource() {
    // Never destroyed due to constant_init!
    return &mpi_res.obj;
}

memory_pool_resource* mpi_pool_resource() {
    // Never destroyed due to constant_init!
    static constant_init<MemoryPoolResource> mpi_pool_res{mpi_resource(), "pluto::mpi_pool_resource", &pluto::memory::mpi_pool};
    return &mpi_pool_res.obj;
}

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto
