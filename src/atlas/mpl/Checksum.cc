/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstring>
#include "atlas/mpl/Checksum.h"

namespace atlas {
namespace mpl {

Checksum::Checksum()
{
  is_setup_ = false;
}

void Checksum::setup( const int part[],
                      const int remote_idx[], const int base,
                      const gidx_t glb_idx[], const gidx_t max_glb_idx,
                      const int parsize )
{
  parsize_ = parsize;
  gather_.setup(part,remote_idx,base,
                glb_idx,max_glb_idx,parsize);
  is_setup_ = true;
}


/////////////////////


Checksum* atlas__Checksum__new () {
  return new Checksum();
}

void atlas__Checksum__delete (Checksum* This) {
  delete This;
}

void atlas__Checksum__setup32 ( Checksum* This,  int part[],
                                int remote_idx[], int base,
                                int glb_idx[], int max_glb_idx,
                                int parsize )
{
#if ATLAS_BITS_GLOBAL==32
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for( int j=0; j<parsize; ++j )
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),max_glb_idx,parsize);
#endif
}

void atlas__Checksum__setup64 ( Checksum* This,  int part[],
                                int remote_idx[], int base,
                                long glb_idx[], long max_glb_idx,
                                int parsize )
{
#if ATLAS_BITS_GLOBAL==64
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for( int j=0; j<parsize; ++j )
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),max_glb_idx,parsize);
#endif
}

void atlas__Checksum__execute_strided_int ( Checksum* This,
                                            int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                            char* checksum )
{
  std::strcpy( checksum, This->execute(lfield,lvar_strides,lvar_extents,lvar_rank).c_str() );
}

void atlas__Checksum__execute_strided_float ( Checksum* This,
                                              float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                              char* checksum )
{
  std::strcpy( checksum, This->execute(lfield,lvar_strides,lvar_extents,lvar_rank).c_str() );
}


void atlas__Checksum__execute_strided_double ( Checksum* This,
                                               double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                               char* checksum )
{
  std::strcpy( checksum, This->execute(lfield,lvar_strides,lvar_extents,lvar_rank).c_str() );
}


//const char* atlas__Checksum__execute_strided_double (Checksum* This,
//                                            double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank)
//{
//  std::string checksum =
//  This->execute(lfield,lvar_strides,lvar_extents,lvar_rank);
//  return checksum.c_str();
//}


/////////////////////

} // namespace mpl
} // namespace mpl
