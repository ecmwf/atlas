#ifndef functionspace_hpp
#define functionspace_hpp

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "Comm.hpp"
#include "Metadata.hpp"

namespace ecmwf {
class Field;
template <typename T> class FieldT;

// TODO:
// Horizontal nodes are always the slowest moving index
// Then variables
// Then levels are fastest moving index
class FunctionSpace
{
public:
  FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& bounds);

  virtual ~FunctionSpace();

  const std::string& name() const { return name_; }

  int index() const { return idx_; }

  Field& field(const std::string& name);

  bool has_field(const std::string& name) { return index_.count(name); }

  template< typename DATA_TYPE>
    FieldT<DATA_TYPE>& field(const std::string& name);

  template< typename DATA_TYPE >
  FieldT<DATA_TYPE>& create_field(const std::string& name, size_t nb_vars);

  void remove_field(const std::string& name);

  const std::vector<int>& bounds() const { return bounds_; }

  void resize( const std::vector<int>& bounds );

  void parallelise(const int proc[], const int glb_idx[], const int master_glb_idx[]);
  void parallelise();

  template< typename DATA_TYPE >
  void halo_exchange( DATA_TYPE field_data[], int field_size )
  {
    int nb_vars = field_size/dof_;
    if( dof_*nb_vars != field_size ) std::cout << "ERROR in FunctionSpace::halo_exchange" << std::endl;
    halo_exchange_.execute( field_data, nb_vars );
  }

  const HaloExchange& halo_exchange() { return halo_exchange_; }

  void set_index(int idx) { idx_ = idx; }

  Metadata& metadata() { return metadata_; }

  template< typename ValueT >
  ValueT& metadata(const std::string name)
  {
    return metadata_.get<ValueT>(name);
  }

  int nb_fields() const { return fields_.size(); }

protected:
  std::string name_;
  int idx_;
  std::vector<int> bounds_;
  std::map< std::string, size_t > index_;
  std::vector< Field* > fields_;
  HaloExchange halo_exchange_;
  int dof_;
  Metadata metadata_;

private:
    // forbid copy constructor by making it private
    FunctionSpace(const FunctionSpace& other);
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
  void ecmwf__FunctionSpace__bounds (FunctionSpace* This, int* &bounds, int &rank);
  Field* ecmwf__FunctionSpace__field (FunctionSpace* This, char* name);
  void ecmwf__FunctionSpace__parallelise (FunctionSpace* This, int proc[], int glb_idx[], int master_glb_idx[]);
  void ecmwf__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size); 
  void ecmwf__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size); 
  void ecmwf__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size); 
  HaloExchange const* ecmwf__FunctionSpace__halo_exchange (FunctionSpace* This); 

}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // functionspace_hpp
