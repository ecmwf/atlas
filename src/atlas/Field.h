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

#ifndef atlas_Field_h
#define atlas_Field_h

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
#include "atlas/Metadata.h"
#include "atlas/Parameters.h"
#include "atlas/State.h"
#include "atlas/util/ArrayView.h"

namespace atlas {

// Forward declaration
template< typename DATA_TYPE > class FieldT;

template<class T>
inline std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
    return eckit::__print_list(s, v);
}

//----------------------------------------------------------------------------------------------------------------------

class Field : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr<Field> Ptr;
  typedef Metadata Parameters;

public: // Static methods

  static Field* create(const eckit::Parametrisation&);
  static Field* create(const ArrayShape&, const eckit::Parametrisation& = Parameters() );
  static Field* create(const Grid&, const eckit::Parametrisation& = Parameters() );

  template<typename T>
  static Field* create(const Array<T>&, const eckit::Parametrisation& = Parameters() );
  static Field* create(const Field&, const eckit::Parametrisation& = Parameters() );

public: // methods

  Field(const eckit::Parametrisation&);

  virtual ~Field();

  Ptr self() { return Ptr(this); }

  template <typename DATA_TYPE> DATA_TYPE* data();
  template <typename DATA_TYPE> DATA_TYPE const* data() const;

  const std::string& data_type() const { return data_type_; }

  const std::string& name() const { return name_; }

  const Metadata& metadata() const { return metadata_; }
  Metadata& metadata() { return metadata_; }

  const std::vector<int>& shapef() const  { return shapef_; }
  const std::vector<size_t>& shape() const { return shape_; }
  const std::vector<size_t>& strides() const { return strides_; }
  size_t stride(int i) const { return strides_[i];}
  size_t shape(int i) const { return shape_[i];}
  size_t nb_vars() const { return nb_vars_; }

  virtual size_t size() const = 0;
  virtual double bytes() const = 0;

  friend std::ostream& operator<<( std::ostream& os, const Field& v);

private: // methods
  
  virtual void print(std::ostream& os) const
  {
      os << "Field[name=" << name()
         << ",datatype=" << data_type_
         << ",size=" << size()
         << ",shape=" << shape_
         << ",strides=" << strides_
         << "]";
  }

private: // members

  virtual void dump( std::ostream& ) const = 0;

  // Allow a State::add() to change the name of the field, only if the field
  // has no name assigned, so it can be used for lookup later.
  friend Field& State::add( Field* );

private: // members

  size_t nb_vars_;

  std::string name_;

  const Grid* grid_;
  FunctionSpace* function_space_;
  Metadata metadata_;


protected: // members
  std::string data_type_;
  std::vector<int> shapef_;
  std::vector<size_t> shape_;
  std::vector<size_t> strides_;

// End of class Field

//------------------------------------------------------------------------------------------------------

/***************************************************************************/
/* Public methods to be removed soon                                       */
/***************************************************************************/

public:


  typedef std::vector< Field::Ptr > Vector;
  enum { UNDEF_VARS = 2147483647 }; // = std::numeric_limits<int>::max() (integer because of fortran)

  Field(const std::string& name, const size_t nb_vars);

  FunctionSpace& function_space() const {
      if( !function_space_ )
          throw eckit::Exception("Field "+name()+" is not associated to any FunctionSpace.");
      return *function_space_;
  }

  virtual void halo_exchange() = 0;

  virtual void allocate(const std::vector<size_t>& shapef)=0;

  const Grid& grid() const {
      if( function_space_ )
      {
          if( function_space_->mesh().has_grid() )
              return function_space_->mesh().grid();
      }
      if( !grid_ )
          throw eckit::Exception("Field "+name()+" is not associated to any Grid.");
      return *grid_;
  }

  const Mesh& mesh() const { ASSERT(function_space_); return function_space_->mesh(); }
  Mesh& mesh() { ASSERT(function_space_); return function_space_->mesh(); }

  void set_function_space(const FunctionSpace& function_space)
  {
    function_space_ = const_cast<FunctionSpace*>(&function_space);
  }

/***************************************************************************/
/* End Public methods to be removed soon                                   */
/***************************************************************************/

};


template< typename DATA_TYPE >
class FieldT : public Field {

public: // methods

  FieldT(const std::string& name, const int nb_vars);

  FieldT(const ArrayShape& shape, const eckit::Parametrisation& = eckit::Properties() );

  virtual ~FieldT();

  virtual size_t size() const { return data_.size(); }

  virtual void allocate(const std::vector<size_t>& shape); // To be removed

  DATA_TYPE* data() { return data_.data(); }
  DATA_TYPE const* data() const { return data_.data(); }

  DATA_TYPE& operator[] (const size_t idx) { return data_[idx]; }

  virtual void halo_exchange(); // To be removed

  virtual void dump(std::ostream& out) const;

  virtual double bytes() const { return sizeof(DATA_TYPE)*size(); }

protected:

  std::vector< DATA_TYPE > data_;

};


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, const int nb_vars) :
  Field(name,nb_vars),
  data_(0)
{
  data_type_ = data_type_to_str<DATA_TYPE>() ;
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::vector<size_t>& shape, const eckit::Parametrisation& params) :
  Field(params),
  data_(0)
{
  data_type_ = data_type_to_str<DATA_TYPE>() ;
  allocate(shape);
}


template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
  data_.clear();
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::allocate(const std::vector<size_t>& shape)
{
  shape_ = shape;
  size_t tot_size(1); for (int i = 0; i < shape_.size(); ++i) tot_size *= shape_[i];
  data_.resize(tot_size);

  shapef_.resize(shape_.size());
  std::reverse_copy( shape_.begin(), shape_.end(), shapef_.begin() );

  strides_.resize(shape_.size());
  strides_[shape_.size()-1] = 1;
  for( int n=shape_.size()-2; n>=0; --n )
  {
    strides_[n] = strides_[n+1]*shape_[n+1];
  }
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::dump(std::ostream& out) const
{
  const FunctionSpace& nodes = function_space();

  ArrayView<DATA_TYPE,1> values( *this );

  ArrayView<DATA_TYPE,2> lonlat( nodes.field("lonlat") );

  ASSERT( values.shape()[0] == lonlat.shape()[0] );

  // Log::info() << "values.shape()[0] " << values.shape()[0] << std::endl;

  for( size_t i = 0; i < lonlat.shape()[0]; ++i )
    out << lonlat(i,LON) << " " << lonlat(i,LAT) << " " << values(i) << std::endl;
}

//----------------------------------------------------------------------------------------------------------------------

#define Parametrisation eckit::Parametrisation
// C wrapper interfaces to C++ routines
extern "C"
{
  Field* atlas__Field__create(Parametrisation* params);
  void atlas__Field__delete (Field* This);
  const char* atlas__Field__name (Field* This);
  const char* atlas__Field__data_type (Field* This);
  int atlas__Field__size (Field* This);
  double atlas__Field__bytes (Field* This);
  int atlas__Field__nb_vars (Field* This);
  void atlas__Field__shapef (Field* This, int* &shape, int &rank);
  void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_shapef, int &rank);
  Metadata* atlas__Field__metadata (Field* This);
  FunctionSpace* atlas__Field__function_space (Field* This);
}
#undef Parametrisation

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
