/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef Checksum_h
#define Checksum_h

#include <vector>

#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/Owned.h"

#include "atlas/mpi/mpi.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Checksum.h"
#include "eckit/utils/Translator.h"

namespace atlas {
namespace mpl {
class Checksum: public eckit::Owned {

public: // types
    typedef eckit::SharedPtr<Checksum> Ptr;

public:
  Checksum();
  Checksum(const std::string& name);
  virtual ~Checksum() {}

public: // methods

  std::string name() const { return name_; }

  /// @brief Setup
  /// @param [in] part         List of partitions
  /// @param [in] remote_idx   List of local indices on remote partitions
  /// @param [in] base         values of remote_idx start at "base"
  /// @param [in] glb_idx      List of global indices
  /// @param [in] max_glb_idx  maximum glb index we want to Checksum.
  ///                          To Checksum everything, set to val > max value in domain
  /// @param [in] parsize      size of given lists
  void setup( const int part[],
              const int remote_idx[], const int base,
              const gidx_t glb_idx[], const gidx_t max_glb_idx,
              const int parsize );

  /// @brief Setup
  /// @param [in] part         List of partitions
  /// @param [in] remote_idx   List of local indices on remote partitions
  /// @param [in] base         values of remote_idx start at "base"
  /// @param [in] glb_idx      List of global indices
  /// @param [in] mask         Mask indices not to include in the communication
  ///                          pattern (0=include,1=exclude)
  /// @param [in] parsize      size of given lists
  void setup( const int part[],
              const int remote_idx[], const int base,
              const gidx_t glb_idx[], const int mask[], const int parsize );

  template <typename DATA_TYPE>
  std::string execute( const DATA_TYPE lfield[],
                       const int lvar_strides[],
                       const int lvar_extents[],
                       const int lvar_rank ) const;

  template <typename DATA_TYPE>
  std::string execute( DATA_TYPE lfield[], const int nb_vars ) const;

  template <typename DATA_TYPE, int LRANK>
  std::string execute( const ArrayView<DATA_TYPE,LRANK>& lfield ) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<int>& varstrides,
                 std::vector<int>& varextents ) const;

private: // data
  std::string name_;
  GatherScatter gather_;
  bool is_setup_;
  int parsize_;
};

template<typename DATA_TYPE>
std::string Checksum::execute( const DATA_TYPE data[],
                               const int var_strides[],
                               const int var_extents[],
                               const int var_rank ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("Checksum was not setup",Here());
  }
  std::vector<checksum_t> local_checksums(parsize_);
  int var_size = var_extents[0]*var_strides[0];

  for( int pp=0; pp<parsize_; ++pp )
  {
    local_checksums[pp] = checksum(data+pp*var_size,var_size);
  }

  std::vector<checksum_t> global_checksums(gather_.glb_dof());
  mpl::Field<checksum_t const> loc(local_checksums.data(),1);
  mpl::Field<checksum_t> glb(global_checksums.data(),1);
  gather_.gather(&loc,&glb,1);

  checksum_t glb_checksum = checksum(global_checksums.data(),global_checksums.size());
  MPI_Bcast(&glb_checksum,1,eckit::mpi::datatype<checksum_t>(),0,eckit::mpi::comm());
  return eckit::Translator<checksum_t,std::string>()(glb_checksum);
}


template <typename DATA_TYPE>
std::string Checksum::execute( DATA_TYPE lfield[], const int nb_vars ) const
{
  int strides[] = {1};
  int extents[] = {nb_vars};
  return execute( lfield, strides, extents, 1 );
}


template<typename DATA_TYPE, int RANK>
void Checksum::var_info( const ArrayView<DATA_TYPE,RANK>& arr,
                         std::vector<int>& varstrides,
                         std::vector<int>& varextents ) const
{
  int rank = std::max(1,RANK-1) ;
  varstrides.resize(rank);
  varextents.resize(rank);
  if( RANK>1 )
  {
    varstrides.assign(arr.strides()+1,arr.strides()+RANK);
    varextents.assign(arr.shape()+1,arr.shape()+RANK);
  }
  else
  {
    varstrides[0] = arr.strides()[0];
    varextents[0] = 1;
  }
}

template <typename DATA_TYPE, int LRANK>
std::string Checksum::execute( const ArrayView<DATA_TYPE,LRANK>& lfield ) const
{
  if( lfield.size() == parsize_ )
  {
    std::vector<int> lvarstrides, lvarextents;
    var_info(lfield, lvarstrides, lvarextents);
    return execute( lfield.data(), lvarstrides.data(), lvarextents.data(), lvarstrides.size() );
  }
  else
  {
    DEBUG_VAR(lfield.size());
    NOTIMP; // Need to implement with parallel ranks > 1
  }
  return std::string("");
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  Checksum* atlas__Checksum__new ();
  void atlas__Checksum__delete (Checksum* This);
  void atlas__Checksum__setup32 (Checksum* This, int part[], int remote_idx[], int base, int glb_idx[], int max_glb_idx, int parsize);
  void atlas__Checksum__setup64 (Checksum* This, int part[], int remote_idx[], int base, long glb_idx[], long max_glb_idx, int parsize);
  void atlas__Checksum__execute_strided_int (Checksum* This, int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, char* checksum);
  void atlas__Checksum__execute_strided_float (Checksum* This, float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, char* checksum);
  void atlas__Checksum__execute_strided_double (Checksum* This, double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, char* checksum);
}
// ------------------------------------------------------------------

} // namespace mpl
} // namespace atlas

#endif // Checksum_h
