
#include <stdexcept>
#include <string>

#define DUMMY_SHOULD_NOT_BE_CALLED(SYMBOL) dummyShouldNotBeCalled( #SYMBOL )
#define DUMMY_FUNCTION(SYMBOL) \
    template <typename... Args> inline \
    dummyError_t dummy##SYMBOL(Args&&... args) { \
        DUMMY_SHOULD_NOT_BE_CALLED( hic##SYMBOL ); \
        return dummyError_t{0}; \
    }
#define DUMMY_VALUE(SYMBOL) \
    constexpr int dummy##SYMBOL = 0;

namespace {

[[noreturn]] void dummyShouldNotBeCalled(const char* symbol) {
    throw std::runtime_error(std::string(symbol)+" is using the dummy backend and should not be called");
}

using dummyError_t  = int;
using dummyEvent_t  = void*;
using dummyStream_t = void*;

struct dummyPointerAttributes {
    int type{0};
    int device{-2};
    void* hostPointer{nullptr};
    void* devicePointer{nullptr};
};
 
DUMMY_FUNCTION(DeviceSynchronize)
DUMMY_FUNCTION(Free)
DUMMY_FUNCTION(FreeAsync)
DUMMY_FUNCTION(GetDeviceCount);
DUMMY_FUNCTION(GetErrorString)
DUMMY_FUNCTION(GetLastError)
DUMMY_FUNCTION(HostGetDevicePointer)
DUMMY_FUNCTION(HostRegister)
DUMMY_FUNCTION(HostUnregister)
DUMMY_FUNCTION(Malloc)
DUMMY_FUNCTION(MallocAsync)
DUMMY_FUNCTION(MallocManaged)
DUMMY_FUNCTION(Memcpy)
DUMMY_FUNCTION(Memcpy2D)
DUMMY_FUNCTION(MemcpyAsync)
DUMMY_FUNCTION(Memcpy2DAsync)
DUMMY_FUNCTION(MemPrefetchAsync)
DUMMY_FUNCTION(StreamCreate)
DUMMY_FUNCTION(StreamDestroy)
DUMMY_FUNCTION(StreamSynchronize)
DUMMY_FUNCTION(PeekAtLastError)
DUMMY_FUNCTION(PointerGetAttributes)

DUMMY_VALUE(CpuDeviceId)
DUMMY_VALUE(HostRegisterMapped)
DUMMY_VALUE(MemoryTypeDevice)
DUMMY_VALUE(MemoryTypeHost)
DUMMY_VALUE(MemoryTypeUnregistered)
DUMMY_VALUE(MemoryTypeManaged)
DUMMY_VALUE(MemcpyDeviceToHost)
DUMMY_VALUE(MemcpyHostToDevice)
DUMMY_VALUE(Success)

}

#undef DUMMY_FUNCTION
#undef DUMMY_VALUE
#undef DUMMY_SHOULD_NOT_BE_CALLED