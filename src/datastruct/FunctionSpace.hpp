#ifndef functionspace_hpp
#define functionspace_hpp

#include <vector>
#include <map>
#include <string>

namespace ecmwf {
class Field;

class FunctionSpace
{
public:
  FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& bounds);

  virtual ~FunctionSpace();
  const std::string& name() const { return name_; }
  Field& field(const std::string& name);

  template< typename DATA_TYPE >
  void create_field(const std::string& name, size_t nb_vars);
  void remove_field(const std::string& name);
  const std::vector<int>& bounds() const { return bounds_; }
protected:
  std::string name_;
  std::vector<int> bounds_;
  std::map< std::string, size_t > index_;
  std::vector< Field* > fields_;
};


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  FunctionSpace* ecmwf__FunctionSpace__new (char* name, char* shape_func, int bounds[], int bounds_size);
  void ecmwf__FunctionSpace__delete (FunctionSpace* This);
  void ecmwf__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars);
  void ecmwf__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars);
  void ecmwf__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars);
  void ecmwf__FunctionSpace__remove_field (FunctionSpace* This, char* name);
  const char* ecmwf__FunctionSpace__name (FunctionSpace* This);
  Field* ecmwf__FunctionSpace__field (FunctionSpace* This, char* name); 
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // functionspace_hpp
