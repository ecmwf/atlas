/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_FunctionSpace_h
#define atlas_FunctionSpace_h

#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/Config.h"

namespace atlas {
namespace next { // Temporary namespace for new design

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpace : public eckit::Owned
{
public:
    FunctionSpace(const std::string& name) : name_(name) {}
    virtual ~FunctionSpace() = 0;
    const std::string& name() const { return name_; }
private:
    std::string name_;
};

inline FunctionSpace::~FunctionSpace() {}

} // namespace next
} // namespace atlas















































//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------


#include <iosfwd>
#include <string>
#include <vector>

#include "eckit/container/DenseMap.h"
#include "eckit/log/Log.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/atlas_config.h"
#include "atlas/mpl/HaloExchange.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/mpl/Checksum.h"
#include "atlas/Metadata.h"
#include "atlas/util/ObjectRegistry.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class Mesh;
class Field;

//------------------------------------------------------------------------------------------------------

enum CreateBehavior { IF_EXISTS_FAIL = 0,    /* when creating, fail if exists */
                      IF_EXISTS_RETURN = 1   /* when creating, return if exists */ };

/// @todo
// Horizontal nodes are always the slowest moving index
// Then variables
// Then levels are fastest moving index
class FunctionSpace : public eckit::Owned, public util::Registered<FunctionSpace> {

public: // types

    typedef eckit::SharedPtr<FunctionSpace> Ptr;

    enum { UNDEF_VARS = 2147483647 }; // = std::numeric_limits<int>::max() (integer because of fortran)

public: // methods

    FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& shape, Mesh& mesh );

    /// TEMPORARY CONSTRUCTOR, JUST FOR EVOLUTIONARY STEP TO NEXT DESIGN
    FunctionSpace(const std::string& name, const std::vector<size_t>& shape );

    virtual ~FunctionSpace();

    const std::string& name() const { return name_; }

    int index() const { return idx_; }

    virtual const Field& field( size_t ) const;
    virtual       Field& field( size_t );
    virtual const Field& field(const std::string& name) const;
    virtual       Field& field(const std::string& name);

    virtual bool has_field(const std::string& name) const { return fields_.has(name); }

    template< typename DATA_TYPE >
    Field& create_field(const std::string& name, size_t nb_vars, CreateBehavior b = IF_EXISTS_FAIL );

    void remove_field(const std::string& name);

    // This is a Fortran view of the shape (i.e. reverse order)
    const std::vector<int>& shapef() const { return shapef_; }

    const std::vector<size_t>& shape() const { return shape_; }
    size_t shape(const size_t i) const { ASSERT(i<shape_.size()); return shape_[i]; }
    void resize( const std::vector<size_t>& shape );


    void parallelise();
    void parallelise(const int proc[], const int remote_idx[], const gidx_t glb_idx[], size_t size);
    void parallelise(FunctionSpace& other_functionspace);

    template< typename DATA_TYPE >
    void halo_exchange( DATA_TYPE field_data[], size_t field_size )
    {
        int nb_vars = field_size/dof();
        if( dof()*nb_vars != field_size )
        {
            eckit::Log::error() << "ERROR in FunctionSpace::halo_exchange" << std::endl;
            eckit::Log::error() << "field_size = " << field_size << std::endl;
            eckit::Log::error() << "dof() = " << dof() << std::endl;
        }
        halo_exchange_->execute( field_data, nb_vars );
    }

    template< typename DATA_TYPE >
    void gather( const DATA_TYPE field_data[], size_t field_size, DATA_TYPE glbfield_data[], size_t glbfield_size )
    {
        int nb_vars = field_size/dof();
        if( dof()*nb_vars != field_size ) eckit::Log::error() << "ERROR in FunctionSpace::gather" << std::endl;
        if( glb_dof_*nb_vars != glbfield_size ) eckit::Log::error() << "ERROR in FunctionSpace::gather" << std::endl;

        mpl::Field<DATA_TYPE const> loc_field(field_data,nb_vars);
        mpl::Field<DATA_TYPE      > glb_field(glbfield_data,nb_vars);

        gather_scatter_->gather( &loc_field, &glb_field, 1 );
    }

    mpl::HaloExchange& halo_exchange() const { return *halo_exchange_; }

    mpl::GatherScatter& gather_scatter() const { return *gather_scatter_; }

    mpl::GatherScatter& fullgather() const { return *fullgather_; }

    mpl::Checksum& checksum() const { return *checksum_; }

    void set_index(size_t idx) { idx_ = idx; }


    const Metadata& metadata() const { return metadata_; }
    Metadata& metadata() { return metadata_; }

    const Mesh& mesh() const { ASSERT(mesh_); return *mesh_; }
    Mesh& mesh() { ASSERT(mesh_); return *mesh_; }

    virtual size_t nb_fields() const { return fields_.size(); }

    size_t dof() const { return dof_; }

    size_t glb_dof() const { return glb_dof_; }

    void print(std::ostream&, bool dump = false) const;

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const FunctionSpace& p) {
        p.print(s);
        return s;
    }

public:
    virtual Field& add( Field* field );

protected:

    size_t dof_;

private:

    size_t idx_;
    size_t glb_dof_;

    std::string name_;

    std::vector<int>    shapef_; // deprecated, use shape which is reverse order
    std::vector<size_t> shape_;

    eckit::DenseMap< std::string, eckit::SharedPtr<Field> > fields_;

    mpl::GatherScatter::Ptr gather_scatter_; // without ghost
    mpl::GatherScatter::Ptr fullgather_; // includes halo
    mpl::HaloExchange::Ptr  halo_exchange_;
    mpl::Checksum::Ptr      checksum_;

    Metadata metadata_;

    Mesh*    mesh_;
};

typedef mpl::HaloExchange HaloExchange_t;
typedef mpl::GatherScatter GatherScatter_t;
typedef mpl::Checksum Checksum_t;

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

extern "C"
{
    Metadata* atlas__FunctionSpace__metadata (FunctionSpace* This);
    int atlas__FunctionSpace__dof (FunctionSpace* This);
    int atlas__FunctionSpace__glb_dof (FunctionSpace* This);
    void atlas__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars);
    void atlas__FunctionSpace__create_field_long (FunctionSpace* This, char* name, int nb_vars);
    void atlas__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars);
    void atlas__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars);
    void atlas__FunctionSpace__remove_field (FunctionSpace* This, char* name);
    int atlas__FunctionSpace__has_field (FunctionSpace* This, char* name);
    const char* atlas__FunctionSpace__name (FunctionSpace* This);
    void atlas__FunctionSpace__shapef (FunctionSpace* This, int* &shape, int &rank);
    Field* atlas__FunctionSpace__field (FunctionSpace* This, char* name);
    void atlas__FunctionSpace__parallelise (FunctionSpace* This);
    void atlas__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size);
    void atlas__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size);
    void atlas__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size);
    void atlas__FunctionSpace__gather_int (FunctionSpace* This, int field_data[], int field_size, int glbfield_data[], int glbfield_size);
    void atlas__FunctionSpace__gather_float (FunctionSpace* This, float field_data[], int field_size, float glbfield_data[], int glbfield_size);
    void atlas__FunctionSpace__gather_double (FunctionSpace* This, double field_data[], int field_size, double glbfield_data[], int glbfield_size);
    HaloExchange_t* atlas__FunctionSpace__halo_exchange (FunctionSpace* This);
    GatherScatter_t* atlas__FunctionSpace__gather (FunctionSpace* This);
    Checksum_t* atlas__FunctionSpace__checksum (FunctionSpace* This);

}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_FunctionSpace_h
