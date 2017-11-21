
#include "atlas/array_fwd.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/StorageView.h"
#include "atlas/library/config.h"

#include "atlas/array.h"


namespace atlas {
namespace array {

namespace {
    template<typename Value, unsigned Rank>
    inline static void check_metadata(const Array& array)
    {
        if(array.rank() != Rank ) {
            std::stringstream err;
            err << "Number of dimensions do not match: template argument " << Rank << " expected to be " << array.rank();
            throw eckit::BadParameter(err.str(), Here());
        }
        if(array.datatype() != array::DataType::create<Value>() ) {
            std::stringstream err;
            err << "Data Type does not match: template argument expected to be " << array.datatype().str();
            throw eckit::BadParameter(err.str(), Here());
        }
    }
}

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode>
make_host_view(const Array& array) {
    return ArrayView<Value, Rank, AccessMode>((const Value*)(array.storage()),array.shape(), array.strides());
}


template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode>
make_device_view(const Array& array) {
    return make_host_view<Value,Rank,AccessMode>(array);
}


template <typename Value, unsigned int Rank, Intent AccessMode>
IndexView<Value, Rank>
make_host_indexview(const Array& array) {
    return IndexView<Value,Rank>( (Value*)(array.storage()),array.shape().data() );
}


template <typename Value>
StorageView<Value>
make_host_storageview(const Array& array) {
    return StorageView<Value>(const_cast<Array&>(array).storage(),array.size(),array.contiguous());
}


template <typename Value>
StorageView<Value>
make_device_storageview(const Array& array) {
    return make_host_storageview<Value>(array);
}


template <typename Value, unsigned int Rank, Intent AccessMode>
IndexView<Value, Rank>
make_indexview(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_indexview<Value,Rank,AccessMode>(array);
}


template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode>
make_view(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_view<Value, Rank, AccessMode>(array);
}


template <typename Value>
StorageView<Value>
make_storageview(const Array& array) {
    return make_host_storageview<Value>(array);
}

// --------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas


//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION(RANK) \
template ArrayView<int,RANK,Intent::ReadOnly> make_view<int,RANK,Intent::ReadOnly>(const Array&);\
template ArrayView<int,RANK,Intent::ReadWrite> make_view<int,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long,RANK,Intent::ReadOnly> make_view<long,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long,RANK,Intent::ReadWrite> make_view<long,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadOnly> make_view<long unsigned,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadWrite> make_view<long unsigned,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<float,RANK,Intent::ReadOnly> make_view<float,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<float,RANK,Intent::ReadWrite> make_view<float,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<double,RANK,Intent::ReadOnly> make_view<double,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<double,RANK,Intent::ReadWrite> make_view<double,RANK,Intent::ReadWrite>(const Array&);\
\
template ArrayView<int,RANK,Intent::ReadOnly> make_host_view<int,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<int,RANK,Intent::ReadWrite> make_host_view<int,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long,RANK,Intent::ReadOnly> make_host_view<long,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long,RANK,Intent::ReadWrite> make_host_view<long,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadOnly> make_host_view<long unsigned,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadWrite> make_host_view<long unsigned,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<float,RANK,Intent::ReadOnly> make_host_view<float,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<float,RANK,Intent::ReadWrite> make_host_view<float,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<double,RANK,Intent::ReadOnly> make_host_view<double,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<double,RANK,Intent::ReadWrite> make_host_view<double,RANK,Intent::ReadWrite>(const Array&);\
\
template ArrayView<int,RANK,Intent::ReadOnly> make_device_view<int,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<int,RANK,Intent::ReadWrite> make_device_view<int,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long,RANK,Intent::ReadOnly> make_device_view<long,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long,RANK,Intent::ReadWrite> make_device_view<long,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadOnly> make_device_view<long unsigned,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<long unsigned,RANK,Intent::ReadWrite> make_device_view<long unsigned,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<float,RANK,Intent::ReadOnly> make_device_view<float,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<float,RANK,Intent::ReadWrite> make_device_view<float,RANK,Intent::ReadWrite>(const Array&);\
template ArrayView<double,RANK,Intent::ReadOnly> make_device_view<double,RANK,Intent::ReadOnly >(const Array&);\
template ArrayView<double,RANK,Intent::ReadWrite> make_device_view<double,RANK,Intent::ReadWrite>(const Array&);\

template IndexView<int,1> make_indexview<int,1,Intent::ReadOnly >(const Array&);\
template IndexView<int,1> make_indexview<int,1,Intent::ReadWrite>(const Array&);\
template IndexView<int,2> make_indexview<int,2,Intent::ReadOnly >(const Array&);\
template IndexView<int,2> make_indexview<int,2,Intent::ReadWrite>(const Array&);\

template IndexView<int,1> make_host_indexview<int,1,Intent::ReadOnly >(const Array&);\
template IndexView<int,1> make_host_indexview<int,1,Intent::ReadWrite>(const Array&);\
template IndexView<int,2> make_host_indexview<int,2,Intent::ReadOnly >(const Array&);\
template IndexView<int,2> make_host_indexview<int,2,Intent::ReadWrite>(const Array&);\

template StorageView<int> make_storageview<int>(const Array&);
template StorageView<long> make_storageview<long>(const Array&);
template StorageView<long unsigned> make_storageview<long unsigned>(const Array&);
template StorageView<float> make_storageview<float>(const Array&);
template StorageView<double> make_storageview<double>(const Array&);

template StorageView<int> make_host_storageview<int>(const Array&);
template StorageView<long> make_host_storageview<long>(const Array&);
template StorageView<long unsigned> make_host_storageview<long unsigned>(const Array&);
template StorageView<float> make_host_storageview<float>(const Array&);
template StorageView<double> make_host_storageview<double>(const Array&);

template StorageView<int> make_device_storageview<int>(const Array&);
template StorageView<long> make_device_storageview<long>(const Array&);
template StorageView<long unsigned> make_device_storageview<long unsigned>(const Array&);
template StorageView<float> make_device_storageview<float>(const Array&);
template StorageView<double> make_device_storageview<double>(const Array&);

// For each Rank in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION(1)
EXPLICIT_TEMPLATE_INSTANTIATION(2)
EXPLICIT_TEMPLATE_INSTANTIATION(3)
EXPLICIT_TEMPLATE_INSTANTIATION(4)
EXPLICIT_TEMPLATE_INSTANTIATION(5)
EXPLICIT_TEMPLATE_INSTANTIATION(6)
EXPLICIT_TEMPLATE_INSTANTIATION(7)
EXPLICIT_TEMPLATE_INSTANTIATION(8)
EXPLICIT_TEMPLATE_INSTANTIATION(9)

#undef EXPLICIT_TEMPLATE_INSTANTIATION
} // namespace array
} // namespace atlas

