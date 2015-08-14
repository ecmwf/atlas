/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_NodesFunctionSpace_h
#define atlas_functionspace_NodesFunctionSpace_h


#include "eckit/memory/SharedPtr.h"
#include "atlas/FunctionSpace.h"
#include "atlas/FieldSet.h"

namespace atlas { class Mesh; }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Halo
{
public:
  Halo(const size_t size) : size_(size) {}
  size_t size() const { return size_; }
private:
  size_t size_;
};

// -------------------------------------------------------------------

class NodesFunctionSpace : public next::FunctionSpace
{
public:

    NodesFunctionSpace(const std::string& name, Mesh& mesh, const Halo& = Halo(0) );

    virtual ~NodesFunctionSpace();

    size_t nb_nodes() const;
    size_t nb_nodes_global() const; // Only on MPI rank 0, will this be different 0

    const Mesh& mesh() const { return mesh_; }
          Mesh& mesh()       { return mesh_; }

    const Nodes& nodes() const { return nodes_; }
          Nodes& nodes()       { return nodes_; }



// -- Local Field creation methods

    /// @brief Create an unnamed scalar field
    template< typename DATATYPE >
    Field* createField() const;

    /// @brief Create an unnamed scalar field
    Field* createField(const std::string& datatype) const;

    /// @brief Create a named scalar field
    template< typename DATATYPE >
    Field* createField(const std::string& name) const;

    /// @brief Create a named scalar field
    Field* createField(const std::string& name, const std::string& datatype) const;

    /// @brief Create an unnamed field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createField(const std::vector<size_t>& variables) const;

    /// @brief Create an unnamed field with specified dimensions for the variables
    Field* createField(const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create an unnamed field based on other field (datatype and dimensioning)
    Field* createField(const Field&) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    Field* createField(const std::string& name, const Field&) const;



// -- Global Field creation methods

    /// @brief Create an unnamed global scalar field
    template< typename DATATYPE >
    Field* createGlobalField() const;

    /// @brief Create an unnamed global scalar field
    Field* createGlobalField(const std::string& datatype) const;

    /// @brief Create a named global scalar field
    template< typename DATATYPE >
    Field* createGlobalField(const std::string& name) const;

    /// @brief Create a named global scalar field
    Field* createGlobalField(const std::string& name, const std::string& datatype) const;

    /// @brief Create an unnamed global field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createGlobalField(const std::vector<size_t>& variables) const;

    /// @brief Create an unnamed global field with specified dimensions for the variables
    Field* createGlobalField(const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create a named global field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create an unnamed global field based on other field (datatype and dimensioning)
    Field* createGlobalField(const Field&) const;

    /// @brief Create a named global field based on other field (datatype and dimensioning)
    Field* createGlobalField(const std::string& name, const Field&) const;



// -- Parallelisation aware methods

    void haloExchange( FieldSet& ) const;
    void haloExchange( Field& ) const;

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;

    template< typename DATATYPE >
    void sum( const Field&, DATATYPE&, size_t& N ) const;

    template< typename DATATYPE >
    void sum( const Field&, std::vector<DATATYPE>&, size_t& N ) const;

    template< typename DATATYPE >
    void orderIndependentSum( const Field&, DATATYPE&, size_t& N ) const;

    template< typename DATATYPE >
    void orderIndependentSum( const Field&, std::vector<DATATYPE>&, size_t& N ) const;

    template< typename DATATYPE >
    void minimum( const Field&, DATATYPE& minimum ) const;

    template< typename DATATYPE >
    void maximum( const Field&, DATATYPE& maximum ) const;

    template< typename DATATYPE >
    void minimum( const Field&, std::vector<DATATYPE>& ) const;

    template< typename DATATYPE >
    void maximum( const Field&, std::vector<DATATYPE>& ) const;

    template< typename DATATYPE >
    void minimumAndLocation( const Field&, DATATYPE& minimum, gidx_t& glb_idx ) const;

    template< typename DATATYPE >
    void maximumAndLocation( const Field&, DATATYPE& maximum, gidx_t& glb_idx ) const;

    template< typename DATATYPE >
    void minimumAndLocation( const Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx ) const;

    template< typename DATATYPE >
    void maximumAndLocation( const Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx ) const;

    template< typename DATATYPE >
    void mean( const Field&, DATATYPE& mean, size_t& N ) const;

    template< typename DATATYPE >
    void mean( const Field&, std::vector<DATATYPE>& mean, size_t& N ) const;

    template< typename DATATYPE >
    void meanAndStandardDeviation( const Field&, DATATYPE& mean, DATATYPE& stddev, size_t& N ) const;

    template< typename DATATYPE >
    void meanAndstandardDeviation( const Field&, std::vector<DATATYPE>& mean, std::vector<DATATYPE>& stddev, size_t& N ) const;

protected: // methods

    std::string halo_name() const;
    std::string gather_scatter_name() const;
    std::string checksum_name() const;

private: // data

    Mesh& mesh_; // non-const because functionspace may modify mesh
    Nodes& nodes_; // non-const because functionspace may modify mesh
    size_t halo_;
    size_t nb_nodes_;
    size_t nb_nodes_global_;

public:
    size_t nb_nodes_global_broadcasted_;
};

// -------------------------------------------------------------------

class NodesColumnFunctionSpace : public NodesFunctionSpace
{
public:

    NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, size_t nb_levels, const Halo& = Halo(0) );

    virtual ~NodesColumnFunctionSpace();

    size_t nb_levels() const;


// -- Local Field creation methods

    /// @brief Create an unnamed scalar field
    template< typename DATATYPE >
    Field* createField() const;

    /// @brief Create an unnamed scalar field
    Field* createField(const std::string& datatype) const;

    /// @brief Create a named scalar field
    template< typename DATATYPE >
    Field* createField(const std::string& name) const;

    /// @brief Create a named scalar field
    Field* createField(const std::string& name, const std::string& datatype) const;

    /// @brief Create an unnamed field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createField(const std::vector<size_t>& variables) const;

    /// @brief Create an unnamed field with specified dimensions for the variables
    Field* createField(const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create an unnamed field based on other field (datatype and dimensioning)
    Field* createField(const Field&) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    Field* createField(const std::string& name, const Field&) const;



// -- Global Field creation methods

    /// @brief Create an unnamed global scalar field
    template< typename DATATYPE >
    Field* createGlobalField() const;

    /// @brief Create an unnamed global scalar field
    Field* createGlobalField(const std::string& datatype) const;

    /// @brief Create a named global scalar field
    template< typename DATATYPE >
    Field* createGlobalField(const std::string& name) const;

    /// @brief Create a named global scalar field
    Field* createGlobalField(const std::string& name, const std::string& datatype) const;

    /// @brief Create an unnamed global field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createGlobalField(const std::vector<size_t>& variables) const;

    /// @brief Create an unnamed global field with specified dimensions for the variables
    Field* createGlobalField(const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create a named global field with specified dimensions for the variables
    template< typename DATATYPE >
    Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

    /// @brief Create an unnamed global field based on other field (datatype and dimensioning)
    Field* createGlobalField(const Field&) const;

    /// @brief Create a named global field based on other field (datatype and dimensioning)
    Field* createGlobalField(const std::string& name, const Field&) const;


// -- Parallelisation aware methods


    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;

    /// @brief Compute sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const Field&, std::vector<DATATYPE>& sum, size_t& N ) const;

    /// @brief Compute sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void sumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute order independent sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute order independent sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const Field&, std::vector<DATATYPE>&, size_t& N ) const;

    /// @brief Compute order independent sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void orderIndependentSumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute minimum of scalar field
    template< typename DATATYPE >
    void minimum( const Field&, DATATYPE& minimum ) const;

    /// @brief Compute maximum of scalar field
    template< typename DATATYPE >
    void maximum( const Field&, DATATYPE& maximum ) const;

    /// @brief Compute minimum of field for each field-variable
    template< typename DATATYPE >
    void minimum( const Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute maximum of field for each field-variable
    template< typename DATATYPE >
    void maximum( const Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute minimum of field for each vertical level separately
    /// @param [out] min    Field of dimension of input without the nodes index
    void minimumPerLevel( const Field&, Field& min ) const;

    /// @brief Compute maximum of field for each vertical level separately
    /// @param [out] max    Field of dimension of input without the nodes index
    void maximumPerLevel( const Field&, Field& max ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, DATATYPE& minimum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, DATATYPE& maximum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute minimum and its location of a field for each vertical level separately
    void minimumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute maximum and its location of a field for each vertical level separately
    void maximumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute mean value of scalar field
    /// @param [out] mean    Mean value
    /// @param [out] N       Number of value used to create the mean
    template< typename DATATYPE >
    void mean( const Field&, DATATYPE& mean, size_t& N ) const;

    /// @brief Compute mean value of field for each field-variable
    /// @param [out] mean    Mean values for each variable
    /// @param [out] N       Number of values used to create the means
    template< typename DATATYPE >
    void mean( const Field&, std::vector<DATATYPE>& mean, size_t& N ) const;

    /// @brief Compute mean values of field for vertical level separately
    /// @param [out] mean    Field of dimension of input without the nodes index
    /// @param [out] N       Number of values used to create the means
    void meanPerLevel( const Field&, Field& mean, size_t& N ) const;

    /// @brief Compute mean value and standard deviation of scalar field
    /// @param [out] mean      Mean value
    /// @param [out] stddev    Standard deviation
    /// @param [out] N         Number of value used to create the mean
    template< typename DATATYPE >
    void meanAndStandardDeviation( const Field&, DATATYPE& mean, DATATYPE& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of scalar field for each field-variable
    /// @param [out] mean      Mean values for each field-variable
    /// @param [out] stddev    Standard deviation for each field-variable
    /// @param [out] N         Number of value used to create the means
    template< typename DATATYPE >
    void meanAndstandardDeviation( const Field&, std::vector<DATATYPE>& mean, std::vector<DATATYPE>& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of field for vertical level separately
    /// @param [out] mean      Field of dimension of input without the nodes index
    /// @param [out] stddev    Field of dimension of input without the nodes index
    /// @param [out] N         Number of values used to create the means
    void meanAndStandardDeviationPerLevel( const Field&, Field& mean, Field& stddev, size_t& N ) const;

private: // data

    size_t nb_levels_;
};

// -------------------------------------------------------------------

template< typename DATATYPE >
Field* NodesFunctionSpace::createField() const
{
    return createField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name) const
{
    return createField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::vector<size_t>& variables) const
{
    return createField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createField(name,variables,DataType::datatype<DATATYPE>());
}


template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField() const
{
    return createGlobalField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name) const
{
    return createGlobalField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::vector<size_t>& variables) const
{
    return createGlobalField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,variables,DataType::datatype<DATATYPE>());
}


template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createField() const
{
    return createField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createField(const std::string& name) const
{
    return createField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createField(const std::vector<size_t>& variables) const
{
    return createField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createField(name,variables,DataType::datatype<DATATYPE>());
}


template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createGlobalField() const
{
    return createGlobalField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name) const
{
  return createGlobalField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createGlobalField(const std::vector<size_t>& variables) const
{
    return createGlobalField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,variables,DataType::datatype<DATATYPE>());
}

#define Char char
extern "C" {
NodesFunctionSpace* atlas__NodesFunctionSpace__new (const char* name, Mesh* mesh, int halo);
void atlas__NodesFunctionSpace__delete (NodesFunctionSpace* This);
int atlas__NodesFunctionSpace__nb_nodes(const NodesFunctionSpace* This);
Mesh* atlas__NodesFunctionSpace__mesh(NodesFunctionSpace* This);
Nodes* atlas__NodesFunctionSpace__nodes(NodesFunctionSpace* This);
Field* atlas__NodesFunctionSpace__create_field (const NodesFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesFunctionSpace__create_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);
Field* atlas__NodesFunctionSpace__create_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template);
Field* atlas__NodesFunctionSpace__create_global_field (const NodesFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesFunctionSpace__create_global_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);
Field* atlas__NodesFunctionSpace__create_global_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template);

void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodesFunctionSpace* This, FieldSet* fieldset);
void atlas__NodesFunctionSpace__halo_exchange_field(const NodesFunctionSpace* This, Field* field);

void atlas__NodesFunctionSpace__gather_fieldset(const NodesFunctionSpace* This, const FieldSet* local, FieldSet* global);
void atlas__NodesFunctionSpace__gather_field(const NodesFunctionSpace* This, const Field* local, Field* global);

void atlas__NodesFunctionSpace__scatter_fieldset(const NodesFunctionSpace* This, const FieldSet* global, FieldSet* local);
void atlas__NodesFunctionSpace__scatter_field(const NodesFunctionSpace* This, const Field* global, Field* local);
void atlas__NodesFunctionSpace__checksum_fieldset(const NodesFunctionSpace* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
void atlas__NodesFunctionSpace__checksum_field(const NodesFunctionSpace* This, const Field* field, Char* &checksum, int &size, int &allocated);

void atlas__NodesFunctionSpace__sum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__sum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__sum_int(const NodesFunctionSpace* This, const Field* field, int &sum, int &N);
void atlas__NodesFunctionSpace__sum_long(const NodesFunctionSpace* This, const Field* field, long &sum, int &N);
void atlas__NodesFunctionSpace__sum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_int(const NodesFunctionSpace* This, const Field* field, int* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_long(const NodesFunctionSpace* This, const Field* field, long* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__oisum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__oisum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__oisum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__oisum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__min_double(const NodesFunctionSpace* This, const Field* field, double &minimum);
void atlas__NodesFunctionSpace__min_float(const NodesFunctionSpace* This, const Field* field, float &minimum);
void atlas__NodesFunctionSpace__min_int(const NodesFunctionSpace* This, const Field* field, int &minimum);
void atlas__NodesFunctionSpace__min_long(const NodesFunctionSpace* This, const Field* field, long &minimum);

void atlas__NodesFunctionSpace__max_double(const NodesFunctionSpace* This, const Field* field, double &maximum);
void atlas__NodesFunctionSpace__max_float(const NodesFunctionSpace* This, const Field* field, float &maximum);
void atlas__NodesFunctionSpace__max_int(const NodesFunctionSpace* This, const Field* field, int &maximum);
void atlas__NodesFunctionSpace__max_long(const NodesFunctionSpace* This, const Field* field, long &maximum);

void atlas__NodesFunctionSpace__min_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, int &size);

void atlas__NodesFunctionSpace__max_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, int &size);

void atlas__NodesFunctionSpace__minloc_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_float(const NodesFunctionSpace* This, const Field* field, float &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_int(const NodesFunctionSpace* This, const Field* field, int &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_long(const NodesFunctionSpace* This, const Field* field, long &minimum, long &glb_idx);

void atlas__NodesFunctionSpace__maxloc_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_float(const NodesFunctionSpace* This, const Field* field, float &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_int(const NodesFunctionSpace* This, const Field* field, int &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_long(const NodesFunctionSpace* This, const Field* field, long &maximum, long &glb_idx);

void atlas__NodesFunctionSpace__minloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__mean_double(const NodesFunctionSpace* This, const Field* field, double &mean, int &N);
void atlas__NodesFunctionSpace__mean_float(const NodesFunctionSpace* This, const Field* field, float &mean, int &N);
void atlas__NodesFunctionSpace__mean_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &N);
void atlas__NodesFunctionSpace__mean_long(const NodesFunctionSpace* This, const Field* field, long &mean, int &N);
void atlas__NodesFunctionSpace__mean_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, int &size, int &N);

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodesFunctionSpace* This, const Field* field, double &mean, double &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodesFunctionSpace* This, const Field* field, float &mean, float &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodesFunctionSpace* This, const Field* field, long &mean, long &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, double* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, float* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, long* &stddev, int &size, int &N);

NodesColumnFunctionSpace* atlas__NodesColumnFunctionSpace__new (const char* name, Mesh* mesh, int nb_levels, int halo);
void atlas__NodesColumnFunctionSpace__delete (NodesColumnFunctionSpace* This);
int atlas__NodesColumnFunctionSpace__nb_levels(const NodesColumnFunctionSpace* This);
Field* atlas__NodesColumnFunctionSpace__create_field (const NodesColumnFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesColumnFunctionSpace__create_field_vars (const NodesColumnFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);
Field* atlas__NodesColumnFunctionSpace__create_field_template (const NodesColumnFunctionSpace* This, const char* name, const Field* field_template);
Field* atlas__NodesColumnFunctionSpace__create_global_field (const NodesColumnFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesColumnFunctionSpace__create_global_field_vars (const NodesColumnFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);
Field* atlas__NodesColumnFunctionSpace__create_global_field_template (const NodesColumnFunctionSpace* This, const char* name, const Field* field_template);
void atlas__NodesColumnFunctionSpace__checksum_fieldset(const NodesColumnFunctionSpace* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
void atlas__NodesColumnFunctionSpace__checksum_field(const NodesColumnFunctionSpace* This, const Field* field, Char* &checksum, int &size, int &allocated);

}
#undef Char

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpace_h
