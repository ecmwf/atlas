#ifndef atlas_GridToolsTraits_h
#define atlas_GridToolsTraits_h

#include "atlas/internals/atlas_defines.h"

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "storage-facility.hpp"
//------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef __CUDACC__
#define BACKEND gridtools::enumtype::Cuda
#else
#define BACKEND gridtools::enumtype::Host
#endif

template <typename Value, unsigned int NDims, bool ReadOnly = false>
using data_view_tt = gridtools::data_view<gridtools::storage_traits<BACKEND>::data_store_t<
                                           Value, gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> >,
                                       ReadOnly>;


}
}
#endif

#endif

