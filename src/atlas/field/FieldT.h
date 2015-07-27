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
#include "eckit/memory/SharedPtr.h"
#include "eckit/config/Parametrisation.h"

#include "atlas/atlas_config.h"
#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Config.h"
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

  FieldT(const ArrayShape& shape, const eckit::Parametrisation& = Config() );

  virtual ~FieldT();

protected:

  virtual size_t size() const { return array_->size(); }

  virtual void halo_exchange(); // To be removed

  virtual double bytes() const { return sizeof(DATA_TYPE)*size(); }

  virtual void dump(std::ostream& os) const
  {
      os << *this << std::endl;
      DATA_TYPE* data = array_->data<DATA_TYPE>();
      for(size_t i = 0; i < array_->size(); ++i)
      {
          os << data[i] << " ";
          if( (i+1)%10 == 0 )
              os << std::endl;
      }
      os << std::endl;
  }

};


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, const int nb_vars) :
  Field(name,nb_vars)
{
  array_.reset( new Array<DATA_TYPE>() );
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::vector<size_t>& shape, const eckit::Parametrisation& params) :
  Field(params)
{
  array_.reset( new Array<DATA_TYPE>() );

  bool fortran(false);
  params.get("fortran",fortran);
  if( fortran ) {
    const std::vector<size_t>& shape_F = shape;
    std::vector<size_t> shape_C (shape_F.size());
    size_t n=shape_F.size();
    for( size_t j=0; j<shape_F.size(); ++j ) shape_C[j] = shape_F[--n];
    resize(shape_C);
  }
  else {
    resize(shape);
  }
}


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
}

//------------------------------------------------------------------------------------------------------

} // namespace field
} // namespace atlas

#endif
