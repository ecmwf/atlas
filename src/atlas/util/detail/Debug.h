#include "eckit/config/Resource.h"
#include "atlas/mesh.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace debug {

inline gidx_t node_global_index(int i=0) {
  static std::vector<gidx_t> g = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_GLOBAL_INDEX", std::vector<gidx_t>{-1} );
  return g[i];
}

inline gidx_t edge_global_index(int i=0) {
  static std::vector<gidx_t> g = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_EDGE_GLOBAL_INDEX", std::vector<gidx_t>{-1} );
  return g[i];
}

inline gidx_t node_uid(int i=0) {
  static std::vector<gidx_t> g = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_UID", std::vector<gidx_t>{-1} );
  return g[i];
}

inline bool is_node_global_index( gidx_t x ) {
  static std::vector<gidx_t> v = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_GLOBAL_INDEX", std::vector<gidx_t>() );
  for( gidx_t g : v ) {
    if ( x == g )
      return true;
  }
  return false;
}

inline bool is_edge_global_index( gidx_t x ) {
  static std::vector<gidx_t> v = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_EDGE_GLOBAL_INDEX", std::vector<gidx_t>() );
  for( gidx_t g : v ) {
    if ( x == g )
      return true;
  }
  return false;
}

inline bool is_node_uid( gidx_t x ) {
  static std::vector<gidx_t> v = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_UID", std::vector<gidx_t>() );
  for( gidx_t g : v ) {
    if ( x == g )
      return true;
  }
  return false;
}

inline int mpi_rank() {
  static int r = eckit::Resource<gidx_t>("$ATLAS_DEBUG_MPI_RANK",-1);
  return r;
}

inline std::string rank_str() {
  static std::string s = "["+std::to_string(parallel::mpi::comm().rank())+"] ";
  return s;
}

}
}
