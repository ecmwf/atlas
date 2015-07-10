/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Tiago Quintino
/// @date Sep 2014

#ifndef atlas_field_FieldT_h
#define atlas_field_FieldT_h


#include <algorithm>
#include <vector>
#include <string>

#include "eckit/types/Types.h"

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/config/Parametrisation.h"

#include "atlas/atlas_config.h"
#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Parametrisation.h"
#include "atlas/Metadata.h"
#include "atlas/Parameters.h"
#include "atlas/State.h"
#include "atlas/Field.h"
#include "atlas/util/DataType.h"


namespace atlas {
namespace field {

//----------------------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class FieldT : public Field {

public: // methods

  FieldT(const std::string& name, const int nb_vars);

  FieldT(const ArrayShape& shape, const eckit::Parametrisation& = Parameters() );

  virtual ~FieldT();

  const DATA_TYPE* data() const { return data_.data(); }
        DATA_TYPE* data()       { return data_.data(); }

  DATA_TYPE& operator[] (const size_t idx) { return data_[idx]; }

  virtual void allocate(const std::vector<size_t>& shape); // To be removed

protected:

  virtual size_t size() const { return data_.size(); }

  virtual void halo_exchange(); // To be removed

  virtual double bytes() const { return sizeof(DATA_TYPE)*size(); }

  virtual void dump(std::ostream& os) const
  {
      os << *this << std::endl;
      for(size_t i = 0; i < data_.size(); ++i)
      {
          os << data_[i] << " ";
          if( (i+1)%10 == 0 )
              os << std::endl;
      }
      os << std::endl;
  }

  typedef std::vector< DATA_TYPE > storage_type;

  storage_type data_;

};


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, const int nb_vars) :
  Field(name,nb_vars),
  data_(0)
{
  data_type_ = DataType::datatype<DATA_TYPE>() ;
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::vector<size_t>& shape, const eckit::Parametrisation& params) :
  Field(params),
  data_(0)
{
  data_type_ = DataType::datatype<DATA_TYPE>() ;

  bool fortran(false);
  params.get("fortran",fortran);
  if( fortran ) {
    std::vector<size_t> shapef(shape.size());
    size_t n=shape.size();
    for( size_t j=0; j<shape.size(); ++j ) shapef[j] = shape[--n];
    allocate(shapef);
  }
  else {
    allocate(shape);
  }
}


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
  data_.clear();
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::allocate(const std::vector<size_t>& shape)
{
    const size_t old_size = data_.size();

    shape_ = shape;
    size_t tot_size = 1;
    for(size_t i = 0; i < shape_.size(); ++i)
        tot_size *= shape_[i];

    if( old_size && tot_size > old_size)
    {
        storage_type new_data(tot_size);
        std::copy(data_.begin(), data_.end(), new_data.begin());
        data_.swap(new_data);
    }
    else
    {
        data_.resize(tot_size);
    }

    shapef_.resize(shape_.size());
    std::reverse_copy( shape_.begin(), shape_.end(), shapef_.begin() );

    strides_.resize(shape_.size());
    strides_[shape_.size()-1] = 1;
    for( int n=shape_.size()-2; n>=0; --n )
    {
        strides_[n] = strides_[n+1]*shape_[n+1];
    }
}

//------------------------------------------------------------------------------------------------------

} // namespace field
} // namespace atlas

#endif
