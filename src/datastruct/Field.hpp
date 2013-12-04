#ifndef field_hpp
#define field_hpp

#include <vector>
#include <string>
#include "Metadata.hpp"
namespace ecmwf {

class FunctionSpace;

class Field {
public:
  enum { NB_VARS = -1 };

  Field(const std::string& name, FunctionSpace& function_space);
  virtual ~Field() {}
  template <typename DATA_TYPE>
    std::vector< DATA_TYPE >& data();
  virtual std::string data_type() const = 0;
  virtual void allocate(const std::vector<int>& bounds)=0;
  const std::string& name() { return name_; }
  Metadata& metadata() { return metadata_; }
  FunctionSpace& function_space() { return function_space_; }
  const std::vector<int>& bounds() const { return bounds_; }
  virtual size_t size() const = 0;
protected:
  std::string name_;
  std::vector<int> bounds_;
  FunctionSpace& function_space_;
  Metadata metadata_;
};

template< typename DATA_TYPE >
class FieldT : public Field
{
public:
  FieldT(const std::string& name, FunctionSpace& function_space);
  virtual ~FieldT();
  virtual std::string data_type() const;
  virtual void allocate(const std::vector<int>& bounds);
  std::vector< DATA_TYPE >& data() { return data_; }
  DATA_TYPE& operator[] (const size_t idx) { return data_[idx]; }
  virtual size_t size() const { return data_.size(); }
protected:
  std::vector< DATA_TYPE > data_;
};

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::FieldT(const std::string& name, FunctionSpace& function_space) :
  Field(name,function_space),
  data_(0)
{
}

template< typename DATA_TYPE >
inline FieldT<DATA_TYPE>::~FieldT()
{
  name_ = "";
  data_.clear();
}

template< typename DATA_TYPE >
inline void FieldT<DATA_TYPE>::allocate(const std::vector<int>& bounds)
{
  bounds_ = bounds;
  size_t tot_size(1); for (int i = 0; i < bounds.size(); ++i) tot_size *= bounds[i];
  data_.resize(tot_size);
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  const char* ecmwf__Field__name (Field* This);
  const char* ecmwf__Field__data_type (Field* This);
  void ecmwf__Field__data_double (Field* This, double* &field_data, int* &field_bounds, int &rank);
  void ecmwf__Field__data_float (Field* This, float* &field_data, int* &field_bounds, int &rank);
  void ecmwf__Field__data_int (Field* This, int* &field_data, int* &field_bounds, int &rank);
  Metadata* ecmwf__Field__metadata (Field* This);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // field_hpp
