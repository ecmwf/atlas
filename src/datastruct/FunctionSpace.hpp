#ifndef functionspace_hpp
#define functionspace_hpp

#include <vector>
#include <map>
#include <string>

namespace ecmwf {
class Field;
class PrismaticField;


class FunctionSpace
{
public:
  FunctionSpace(const std::string& name, const std::string& shape_func, size_t nb_nodes);
  FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& bounds);

  virtual ~FunctionSpace();
  const std::string& name() const { return name_; }
  Field& field(const std::string& name);
  virtual void create_field(const std::string& name, size_t nb_cols);
  void remove_field(const std::string& name);
  void print_field(const std::string& name);
  const std::vector<size_t>& bounds() const { return bounds_; }
protected:
  std::string name_;
  std::vector<size_t> bounds_;
  std::map< std::string, size_t > index_;
  std::vector< Field* > fields_;
};

class PrismaticFunctionSpace: public FunctionSpace
{
public:
  PrismaticFunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& bounds);
  virtual ~PrismaticFunctionSpace();
  virtual void create_field(const std::string& name, size_t nb_cols);
  PrismaticField& field(const std::string& name);
  size_t nb_levels() const { return bounds_[0]; }
  size_t nb_nodes()  const { return bounds_[1]; }

private:
  size_t nb_levels_;
  size_t nb_nodes_;
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  PrismaticFunctionSpace* ecmwf__PrismaticFunctionSpace__new (char* name, char* shape_func, int nb_levels, int nb_nodes);
  FunctionSpace* ecmwf__FunctionSpace__new (char* name, char* shape_func, int nb_nodes);
  void ecmwf__FunctionSpace__create_field (FunctionSpace* This, char* name, int nb_cols);
  void ecmwf__FunctionSpace__remove_field (FunctionSpace* This, char* name);
  const char* ecmwf__FunctionSpace__name (FunctionSpace* This);
  Field* ecmwf__FunctionSpace__field (FunctionSpace* This, char* name); 
  void ecmwf__FunctionSpace__field_data (FunctionSpace* This, char* name, double* &field_data, int &field_size); 
  void ecmwf__FunctionSpace__print_field (FunctionSpace* This, char* name);
  void ecmwf__FunctionSpace__delete (FunctionSpace* This);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // functionspace_hpp
