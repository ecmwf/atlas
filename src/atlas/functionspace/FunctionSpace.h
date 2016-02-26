/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/util/Config.h"

namespace atlas {
namespace functionspace {

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpace : public eckit::Owned
{
public:
    FunctionSpace() {}
    virtual ~FunctionSpace() = 0;
    virtual std::string name() const = 0;
    eckit::SharedPtr<FunctionSpace const> shared_from_this() const;
    eckit::SharedPtr<FunctionSpace> shared_from_this();
    eckit::SharedPtr<FunctionSpace> ptr();
    eckit::SharedPtr<FunctionSpace const> ptr() const;
    eckit::SharedPtr<FunctionSpace const> cptr() const;
};

inline FunctionSpace::~FunctionSpace() {}

} // namespace functionspace
} // namespace atlas


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

// The *old* FunctionSpace class is still in use for the element-types (quads, triags, edges)
// The nodes are now in the class Nodes, which for now still inherits from following.

#include <iosfwd>
#include <string>
#include <vector>

#include "eckit/container/DenseMap.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/atlas_config.h"
#include "atlas/util/runtime/Log.h"
#include "atlas/util/parallel/mpl/HaloExchange.h"
#include "atlas/util/parallel/mpl/GatherScatter.h"
#include "atlas/util/parallel/mpl/Checksum.h"
#include "atlas/util/Metadata.h"
#include "atlas/internals/ObjectRegistry.h"

//------------------------------------------------------------------------------------------------------
namespace atlas { namespace mesh { class Mesh; } }
namespace atlas { namespace field { class Field; } }

namespace atlas {
namespace functionspace {


//------------------------------------------------------------------------------------------------------

namespace deprecated {
#if !DEPRECATE_OLD_FUNCTIONSPACE

enum CreateBehavior { IF_EXISTS_FAIL = 0,    /* when creating, fail if exists */
                      IF_EXISTS_RETURN = 1   /* when creating, return if exists */ };

/// @todo
// Horizontal nodes are always the slowest moving index
// Then variables
// Then levels are fastest moving index
class FunctionSpace : public eckit::Owned, public internals::Registered<FunctionSpace> {

public: // types

    typedef eckit::SharedPtr<FunctionSpace> Ptr;

    enum { UNDEF_VARS = 2147483647 }; // = std::numeric_limits<int>::max() (integer because of fortran)

public: // methods

    FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& shape, mesh::Mesh& mesh );

    /// TEMPORARY CONSTRUCTOR, JUST FOR EVOLUTIONARY STEP TO NEXT DESIGN
    FunctionSpace(const std::string& name, const std::vector<size_t>& shape );

    virtual ~FunctionSpace();

    const std::string& name() const { return name_; }

    int index() const { return idx_; }

    virtual const field::Field& field( size_t ) const;
    virtual       field::Field& field( size_t );
    virtual const field::Field& field(const std::string& name) const;
    virtual       field::Field& field(const std::string& name);

    virtual bool has_field(const std::string& name) const { return fields_.has(name); }

    template< typename DATA_TYPE >
    field::Field& create_field(const std::string& name, size_t nb_vars, CreateBehavior b = IF_EXISTS_FAIL );

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
            Log::error() << "ERROR in FunctionSpace::halo_exchange" << std::endl;
            Log::error() << "field_size = " << field_size << std::endl;
            Log::error() << "dof() = " << dof() << std::endl;
        }
        halo_exchange_->execute( field_data, nb_vars );
    }

    template< typename DATA_TYPE >
    void gather( const DATA_TYPE field_data[], size_t field_size, DATA_TYPE glbfield_data[], size_t glbfield_size )
    {
        int nb_vars = field_size/dof();
        if( dof()*nb_vars != field_size ) Log::error() << "ERROR in FunctionSpace::gather" << std::endl;
        if( glb_dof_*nb_vars != glbfield_size ) Log::error() << "ERROR in FunctionSpace::gather" << std::endl;

        util::parallel::mpl::field::Field<DATA_TYPE const> loc_field(field_data,nb_vars);
        util::parallel::mpl::field::Field<DATA_TYPE      > glb_field(glbfield_data,nb_vars);

        gather_scatter_->gather( &loc_field, &glb_field, 1 );
    }

    util::parallel::mpl::HaloExchange& halo_exchange() const { return *halo_exchange_; }

    util::parallel::mpl::GatherScatter& gather_scatter() const { return *gather_scatter_; }

    util::parallel::mpl::GatherScatter& fullgather() const { return *fullgather_; }

    util::parallel::mpl::Checksum& checksum() const { return *checksum_; }

    void set_index(size_t idx) { idx_ = idx; }


    const util::Metadata& metadata() const { return metadata_; }
    util::Metadata& metadata() { return metadata_; }

    const mesh::Mesh& mesh() const { ASSERT(mesh_); return *mesh_; }
    mesh::Mesh& mesh() { ASSERT(mesh_); return *mesh_; }

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
    virtual field::Field& add( field::Field* field );

protected:

    size_t dof_;

private:

    size_t idx_;
    size_t glb_dof_;

    std::string name_;

    std::vector<int>    shapef_; // deprecated, use shape which is reverse order
    std::vector<size_t> shape_;

    eckit::DenseMap< std::string, eckit::SharedPtr<field::Field> > fields_;

    util::parallel::mpl::GatherScatter::Ptr gather_scatter_; // without ghost
    util::parallel::mpl::GatherScatter::Ptr fullgather_; // includes halo
    util::parallel::util::parallel::mpl::HaloExchange::Ptr  halo_exchange_;
    util::parallel::mpl::Checksum::Ptr      checksum_;

    util::Metadata metadata_;

    mesh::Mesh*    mesh_;
};
#endif
} // namespace deprecated

typedef util::parallel::mpl::HaloExchange HaloExchange_t;
typedef util::parallel::mpl::GatherScatter GatherScatter_t;
typedef util::parallel::mpl::Checksum Checksum_t;

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define util_Metadata util::Metadata
#define field_Field field::Field

#if !DEPRECATE_OLD_FUNCTIONSPACE
#define DeprecatedFunctionSpace deprecated::FunctionSpace
#else
#define DeprecatedFunctionSpace void
#endif
extern "C"
{
    void atlas__FunctionSpace__delete (FunctionSpace* This);
    const char* atlas__FunctionSpace__name (FunctionSpace* This);
    util_Metadata* atlas__deprecated__FunctionSpace__metadata (DeprecatedFunctionSpace* This);
    int atlas__deprecated__FunctionSpace__dof (DeprecatedFunctionSpace* This);
    int atlas__deprecated__FunctionSpace__glb_dof (DeprecatedFunctionSpace* This);
    void atlas__deprecated__FunctionSpace__create_field_int (DeprecatedFunctionSpace* This, char* name, int nb_vars);
    void atlas__deprecated__FunctionSpace__create_field_long (DeprecatedFunctionSpace* This, char* name, int nb_vars);
    void atlas__deprecated__FunctionSpace__create_field_float (DeprecatedFunctionSpace* This, char* name, int nb_vars);
    void atlas__deprecated__FunctionSpace__create_field_double (DeprecatedFunctionSpace* This, char* name, int nb_vars);
    void atlas__deprecated__FunctionSpace__remove_field (DeprecatedFunctionSpace* This, char* name);
    int atlas__deprecated__FunctionSpace__has_field (DeprecatedFunctionSpace* This, char* name);
    const char* atlas__deprecated__FunctionSpace__name (DeprecatedFunctionSpace* This);
    void atlas__deprecated__FunctionSpace__shapef (DeprecatedFunctionSpace* This, int* &shape, int &rank);
    field_Field* atlas__deprecated__FunctionSpace__field (DeprecatedFunctionSpace* This, char* name);
    void atlas__deprecated__FunctionSpace__parallelise (DeprecatedFunctionSpace* This);
    void atlas__deprecated__FunctionSpace__halo_exchange_int (DeprecatedFunctionSpace* This, int field_data[], int field_size);
    void atlas__deprecated__FunctionSpace__halo_exchange_float (DeprecatedFunctionSpace* This, float field_data[], int field_size);
    void atlas__deprecated__FunctionSpace__halo_exchange_double (DeprecatedFunctionSpace* This, double field_data[], int field_size);
    void atlas__deprecated__FunctionSpace__gather_int (DeprecatedFunctionSpace* This, int field_data[], int field_size, int glbfield_data[], int glbfield_size);
    void atlas__deprecated__FunctionSpace__gather_float (DeprecatedFunctionSpace* This, float field_data[], int field_size, float glbfield_data[], int glbfield_size);
    void atlas__deprecated__FunctionSpace__gather_double (DeprecatedFunctionSpace* This, double field_data[], int field_size, double glbfield_data[], int glbfield_size);
    HaloExchange_t* atlas__deprecated__FunctionSpace__halo_exchange (DeprecatedFunctionSpace* This);
    GatherScatter_t* atlas__deprecated__FunctionSpace__gather (DeprecatedFunctionSpace* This);
    Checksum_t* atlas__deprecated__FunctionSpace__checksum (DeprecatedFunctionSpace* This);

}
#undef field_Field
#undef util_Metadata
#undef DeprecatedFunctionSpace
//------------------------------------------------------------------------------------------------------

} //namespace functionspace
} // namespace atlas

#endif // atlas_FunctionSpace_h
