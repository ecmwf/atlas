#ifndef field_hpp
#define field_hpp

#include <vector>
#include <string>
#include "Metadata.hpp"
namespace ecmwf {

class FunctionSpace;
class PrismaticFunctionSpace;

class Field
{
public:
  Field(const std::string& name, const size_t nb_cols, FunctionSpace& function_space);
  virtual ~Field();
  size_t size() const { return data_.size(); }
  std::vector<double>& data() { return data_; }
  double& operator[] (const size_t idx) { return data_[idx]; }
  const std::string& name() { return name_; }
  Metadata& metadata() { return metadata_; }
  virtual void allocate();
  FunctionSpace& function_space() { return function_space_; }
  const std::vector<int>& bounds() const { return bounds_; }
protected:
  std::string name_;
  std::vector<int> bounds_;
  std::vector<double> data_;
  int nb_cols_;
  FunctionSpace& function_space_;
  Metadata metadata_;
};

class PrismaticField : public Field
{
public:
  PrismaticField(const std::string& name, const size_t nb_cols, FunctionSpace& function_space);
  ~PrismaticField();
  double& operator() (size_t lev, size_t node, size_t col) { return data_[lev + nb_levels_*(node + nb_nodes_*col)]; }
  PrismaticFunctionSpace& function_space();
private:
  int nb_nodes_;
  int nb_levels_;
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  const char* ecmwf__Field__name (Field* This);
  void ecmwf__Field__data (Field* This, double* &field_data, int* &field_bounds, int &rank);
  Metadata* ecmwf__Field__metadata (Field* This);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // field_hpp
