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

#ifndef atlas_mesh_Field_hpp
#define atlas_mesh_Field_hpp

#include <algorithm>
#include <vector>
#include <string>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/ScopedPtr.h"

#include "atlas/atlas_config.h"

#ifdef ECKIT_HAVE_GRIB
  #include "eckit/grib/GribHandle.h"
#endif

#include "atlas/grid/Grid.h"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/util/ArrayView.hpp"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

class Field : public eckit::Owned {

public: // types

	typedef eckit::SharedPtr<Field> Ptr;
	typedef std::vector< Field::Ptr > Vector;
	typedef grid::Grid Grid;

#ifdef ECKIT_HAVE_GRIB
	typedef eckit::grib::GribHandle Grib;
#endif

public: // methods

  enum { UNDEF_VARS = -1 };

  Field(const std::string& name, const int nb_vars, FunctionSpace& function_space);

  virtual ~Field();

  template <typename DATA_TYPE> DATA_TYPE* data();
  template <typename DATA_TYPE> DATA_TYPE const* data() const;

  const std::string& data_type() const { return data_type_; }

  virtual void allocate(const std::vector<int>& bounds)=0;
  const std::string& name() const { return name_; }

  const Grid& grid() const { return mesh().grid(); }
  Grid& grid() { return mesh().grid(); }

  const Mesh& mesh() const { return function_space_.mesh(); }
  Mesh& mesh() { return function_space_.mesh(); }

  const Metadata& metadata() const { return metadata_; }
  Metadata& metadata() { return metadata_; }

  FunctionSpace& function_space() const { return function_space_; }

  const std::vector<int>& boundsf() const  { return bounds_; }
  const std::vector<int>& extents() const { return extents_; }
  const std::vector<int>& strides() const { return strides_; }
  int stride(int i) const { return strides_[i];}
  int extent(int i) const { return extents_[i];}
  int nb_vars() const { return nb_vars_; }

  virtual size_t size() const = 0;
  virtual void halo_exchange() = 0;

#ifdef ECKIT_HAVE_GRIB
  void grib( Grib* g );   ///< @todo this is to be removed
  Grib& grib() const;     ///< @todo this is to be removed
#endif

  friend std::ostream& operator<<( std::ostream& os, const Field& v);

private: // members

	virtual void print( std::ostream& ) const = 0;

protected: // members

  std::string name_;
  std::string data_type_;
  std::vector<int> bounds_;
  std::vector<int> extents_;
  std::vector<int> strides_;

  FunctionSpace& function_space_;
  Metadata metadata_;

  int nb_vars_;

#ifdef ECKIT_HAVE_GRIB
  eckit::ScopedPtr<Grib> grib_; ///< @todo this is to be removed
#endif

};

//------------------------------------------------------------------------------------------------------


template< typename DATA_TYPE > inline std::string data_type_to_str();
template<> inline std::string data_type_to_str<int>()    { return "int32";  }
template<> inline std::string data_type_to_str<float>()  { return "real32"; }
template<> inline std::string data_type_to_str<double>() { return "real64"; }

template< typename DATA_TYPE >
class FieldT : public Field {

public: // methods

    FieldT(const std::string& name, const int nb_vars, FunctionSpace& function_space);

    virtual ~FieldT();

    virtual size_t size() const { return data_.size(); }

    virtual void allocate(const std::vector<int>& bounds);

    DATA_TYPE* data() { return data_.data(); }
    DATA_TYPE const* data() const { return data_.data(); }

    DATA_TYPE& operator[] (const size_t idx) { return data_[idx]; }

    virtual void halo_exchange();

	virtual void print(std::ostream& out) const;

protected:

    std::vector< DATA_TYPE > data_;

};

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, const int nb_vars, FunctionSpace& function_space) :
  Field(name,nb_vars,function_space),
  data_(0)
{
  data_type_ = data_type_to_str<DATA_TYPE>() ;
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
    data_.clear();
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::allocate(const std::vector<int>& extents)
{
  extents_ = extents;
  size_t tot_size(1); for (int i = 0; i < extents_.size(); ++i) tot_size *= extents_[i];
  data_.resize(tot_size);

  bounds_.resize(extents_.size());
  std::reverse_copy( extents_.begin(), extents_.end(), bounds_.begin() );

  strides_.resize(extents_.size());
  strides_[extents_.size()-1] = 1;
  for( int n=extents_.size()-2; n>=0; --n )
  {
    strides_[n] = strides_[n+1]*extents_[n+1];
  }
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::print(std::ostream& out) const
{
	const FunctionSpace& nodes = function_space();

	ArrayView<DATA_TYPE,1> values( *this );

//    ArrayView<DATA_TYPE,2> coords( nodes.field("coordinates") );
	ArrayView<DATA_TYPE,2> latlon( nodes.field("latlon") );

	ASSERT( values.extents()[0] == latlon.extents()[0] );

//    Log::info() << "values.extents()[0] " << values.extents()[0] << std::endl;

	for( size_t i = 0; i < latlon.extents()[0]; ++i )
		out << latlon(i,LAT) << " " << latlon(i,LON) << " " << values(i) << std::endl;
}

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" 
{
  const char* atlas__Field__name (Field* This);
  const char* atlas__Field__data_type (Field* This);
  int atlas__Field__nb_vars (Field* This);
  void atlas__Field__data_boundsf_double (Field* This, double* &field_data, int* &field_bounds, int &rank);
  void atlas__Field__data_boundsf_float (Field* This, float* &field_data, int* &field_bounds, int &rank);
  void atlas__Field__data_boundsf_int (Field* This, int* &field_data, int* &field_bounds, int &rank);
  Metadata* atlas__Field__metadata (Field* This);
  FunctionSpace* atlas__Field__function_space (Field* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
