// (C) Copyright 1996-2014 ECMWF.

#ifndef atlas_Field_hpp
#define atlas_Field_hpp

#include <algorithm>
#include <vector>
#include <string>

#include "atlas/Metadata.hpp"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class FunctionSpace;

//------------------------------------------------------------------------------------------------------

class Field {

public: // methods

  enum { UNDEF_VARS = -1 };

  Field(const std::string& name, const int nb_vars, FunctionSpace& function_space);

  virtual ~Field();

  template <typename DATA_TYPE> std::vector< DATA_TYPE >& data();
  template <typename DATA_TYPE> const std::vector< DATA_TYPE >& data() const;

  virtual std::string data_type() const = 0;

  virtual void allocate(const std::vector<int>& bounds)=0;
  const std::string& name() const { return name_; }

  Metadata& metadata() { return metadata_; }

  FunctionSpace& function_space() { return function_space_; }

  const std::vector<int>& bounds() const  { return bounds_; }
  const std::vector<int>& extents() const { return extents_; }
  const std::vector<int>& strides() const { return strides_; }
  int stride(int i) const { return strides_[i];}
  int extent(int i) const { return extents_[i];}
  int nb_vars() const { return nb_vars_; }

  virtual size_t size() const = 0;
  virtual void halo_exchange() = 0;

protected: // members

  std::string name_;
  std::vector<int> bounds_;
  std::vector<int> extents_;
  std::vector<int> strides_;
  FunctionSpace& function_space_;
  Metadata metadata_;
  int nb_vars_;

private: // copy not allowed

    Field(const Field&);
    Field& operator=(const Field&);
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class FieldT : public Field {

public: // methods

    FieldT(const std::string& name, const int nb_vars, FunctionSpace& function_space);

    virtual ~FieldT();

    virtual size_t size() const { return data_.size(); }

    virtual std::string data_type() const;

    virtual void allocate(const std::vector<int>& bounds);

    std::vector< DATA_TYPE >& data() { return data_; }
    const std::vector< DATA_TYPE >& data() const { return data_; }

    DATA_TYPE& operator[] (const size_t idx) { return data_[idx]; }

    DATA_TYPE& operator() (int i, int j)
    {
        return *(data_.data()+ i + j*nb_vars_);
    }

    DATA_TYPE& operator() (int i)
    {
        return data_[i];
    }

    DATA_TYPE* slice( int j )
    {
        return data_.data() + j*nb_vars_;
    }

    virtual void halo_exchange();

protected:

    std::vector< DATA_TYPE > data_;

};

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, const int nb_vars, FunctionSpace& function_space) :
  Field(name,nb_vars,function_space),
  data_(0)
{
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
    data_.clear();
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::allocate(const std::vector<int>& bounds)
{
  bounds_ = bounds;
  size_t tot_size(1); for (int i = 0; i < bounds_.size(); ++i) tot_size *= bounds_[i];
  data_.resize(tot_size);
  
  extents_.resize(bounds_.size());
  std::reverse_copy( bounds_.begin(), bounds_.end(), extents_.begin() );
  strides_.resize(extents_.size());
  strides_[extents_.size()-1] = 1;
  for( int n=extents_.size()-2; n>=0; --n )
  {
    strides_[n] = strides_[n+1]*extents_[n+1];
  }
}

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" 
{
  const char* atlas__Field__name (Field* This);
  const char* atlas__Field__data_type (Field* This);
  int atlas__Field__nb_vars (Field* This);
  void atlas__Field__data_double (Field* This, double* &field_data, int* &field_bounds, int &rank);
  void atlas__Field__data_float (Field* This, float* &field_data, int* &field_bounds, int &rank);
  void atlas__Field__data_int (Field* This, int* &field_data, int* &field_bounds, int &rank);
  Metadata* atlas__Field__metadata (Field* This);
  FunctionSpace* atlas__Field__function_space (Field* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
