#ifndef functionspace_hpp
#define functionspace_hpp

#include <vector>
#include <map>
#include <string>
#include "Comm.hpp"
#include <iostream>

namespace ecmwf {
class Field;

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
  Field& field(const std::string& name);

  template< typename DATA_TYPE >
  void create_field(const std::string& name, size_t nb_vars);
  void remove_field(const std::string& name);
  const std::vector<int>& bounds() const { return bounds_; }

  void parallelise(const int proc[], const int glb_idx[]);

  template< typename DATA_TYPE >
  void halo_exchange( DATA_TYPE field_data[], int field_size )
  {
    int nb_vars = field_size/dof_;
    if( dof_*nb_vars != field_size ) std::cout << "ERROR in FunctionSpace::halo_exchange" << std::endl;
    halo_exchange_.execute( field_data, nb_vars );
  }

  const HaloExchange& halo_exchange() { return halo_exchange_; }

protected:
  std::string name_;
  std::vector<int> bounds_;
  std::map< std::string, size_t > index_;
  std::vector< Field* > fields_;
  int par_extent_; ///! The dimension which is parallelised (e.g. the horizontal node-idx)
  HaloExchange halo_exchange_;
  int dof_;
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
  void ecmwf__FunctionSpace__parallelise (FunctionSpace* This, int proc[], int glb_idx[]);
  void ecmwf__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size); 
  void ecmwf__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size); 
  void ecmwf__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size); 
  HaloExchange const* ecmwf__FunctionSpace__halo_exchange (FunctionSpace* This); 

}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // functionspace_hpp
