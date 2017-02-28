/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/functionspace/StructuredColumns.h"

#include "eckit/utils/MD5.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/internals/Checksum.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/mpi/mpi.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif

namespace atlas {
namespace functionspace {


namespace {
void set_field_metadata(const eckit::Parametrisation& config, field::Field& field)
{
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      field.metadata().set("owner",owner);
    }
  }
  field.metadata().set("global",global);
}
}

size_t StructuredColumns::config_size(const eckit::Parametrisation& config) const
{
  size_t size = npts();
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      size = (eckit::mpi::rank() == owner ? grid_.npts() : 0);
    }
  }
  return size;
}


// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
StructuredColumns::StructuredColumns(const grid::Grid& grid) :
  FunctionSpace(),
  grid_(grid)
{
    if ( not grid_ )
    {
      throw eckit::BadCast("Grid is not a grid::Structured type. "
                           "Cannot partition using IFS trans", Here());
    }

#ifdef ATLAS_HAVE_TRANS
    trans_ = new trans::Trans(grid_);

    npts_ = trans_->ngptot();

    // Number of latitude bands
    int n_regions_NS = trans_->n_regions_NS();

    // Number of partitions per latitude band
    array::ArrayView<int,1> n_regions = trans_->n_regions();

    // First latitude of latitude band
    array::ArrayView<int,1> nfrstlat = trans_->nfrstlat();

    // First latitude of latitude band
    array::ArrayView<int,1> nlstlat = trans_->nlstlat();

    // Index of latitude partition (note that if a partition
    // has two regions on a latitude - the index increases
    // by one (2 numbers)
    array::ArrayView<int,1> nptrfrstlat = trans_->nptrfrstlat();

    // Starting longitudinal point per given latitude (ja)
    // Note that it is associated to nptrfrstlat
    array::ArrayView<int,2> nsta = trans_->nsta();

    // Number of longitudinal points per given latitude (ja)
    // Note that it is associated to nptrfrstlat
    array::ArrayView<int,2> nonl = trans_->nonl();

    size_t proc(0);
    // Loop over number of latitude bands (ja)
    for (int ja = 0; ja < n_regions_NS; ++ja)
    {
        // Loop over number of longitude bands (jb)
        for (int jb = 0; jb < n_regions[ja]; ++jb)
        {
            if (proc == eckit::mpi::rank())
            {
                nlat_ = nlstlat[ja] - nfrstlat[ja] + 1;
                nlon_.resize(nlat_);
                first_lon_.resize(nlat_);
                first_lat_ = nfrstlat[ja]-1;

                // Loop over latitude points of lat band (ja) and lon band (jb)
                size_t ilat = 0;
                for (int jglat = first_lat_; jglat < nlstlat[ja]; ++jglat)
                {
                    size_t igl = nptrfrstlat[ja] + jglat - nfrstlat[ja];
                    nlon_[ilat] = nonl(jb,igl);
                    first_lon_[ilat] = nsta(jb,igl);
                    ilat++;
                }
                goto exit_outer_loop;
            }
            ++proc;
        }
    }
    exit_outer_loop: ;


#if 0
    // COMPILED OUT
    {
        int localLatID;
        int localLonID;
        proc = 0;

        // Loop over number of latitude bands (ja)
        for (size_t ja = 0; ja < n_regions_NS; ++ja)
        {
            // Loop over number of longitude bands (jb)
            for (size_t jb = 0; jb < n_regions[ja]; ++jb)
            {
                if (proc == eckit::mpi::rank())
                {
                    // Loop over latitude points of lat band (ja) and lon band (jb)
                    for (int jglat = nfrstlat[ja]-1; jglat < nlstlat[ja]; ++jglat)
                    {
                        int globalLatID = localLatID + nfrstlat[ja];
                        size_t igl = nptrfrstlat[ja] + jglat - nfrstlat[ja];

                        // Loop over longitude points of given latitude point
                        // of lat band (ja) and lon band (jb) and
                        for (int jglon = nsta(jb,igl)-1;
                             jglon < nsta(jb,igl)+nonl(jb,igl)-1; ++jglon)
                        {
                            int globalLonID = nsta(jb,igl) + localLonID;
                        }
                    }
                }
                ++proc;
            }
        }
    }
#endif
#else
    npts_ = grid_.npts();
    nlat_ = grid_.ny();
    nlon_.resize(nlat_);
    for( size_t jlat=0; jlat<nlat_; ++jlat )
      nlon_[jlat] = grid_.nx(jlat);
    first_lat_ = 0;
    first_lon_.resize(nlat_,0);
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------
StructuredColumns::~StructuredColumns()
{
#ifdef ATLAS_HAVE_TRANS
    delete trans_;
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
field::Field* StructuredColumns::createField(const std::string& name, array::DataType datatype, const eckit::Parametrisation& options ) const
{
#ifdef ATLAS_HAVE_TRANS
    size_t npts = config_size(options);
    field::Field* field = field::Field::create(name, datatype, array::make_shape(npts));
    field->set_functionspace(*this);
    set_field_metadata(options,*field);
    return field;
#else
    if( eckit::mpi::size() > 1 )
    {
        throw eckit::NotImplemented(
          "StructuredColumns::createField currently relies"
          " on ATLAS_HAVE_TRANS for parallel fields", Here());
    }
    field::Field* field = field::Field::create(name, datatype, array::make_shape(grid_.npts()) );
    field->set_functionspace(*this);
    set_field_metadata(options,*field);
    return field;
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create Field with vertical levels
// ----------------------------------------------------------------------------
field::Field* StructuredColumns::createField(
    const std::string& name, array::DataType datatype,
    size_t levels, const eckit::Parametrisation& options) const
{
#ifdef ATLAS_HAVE_TRANS
    size_t npts = config_size(options);
    field::Field* field = field::Field::create<double>(
                    name, array::make_shape(npts, levels));

    field->set_functionspace(*this);
    field->set_levels(levels);
    set_field_metadata(options,*field);
    return field;
#else
    if( eckit::mpi::size() > 1 )
    {
        throw eckit::NotImplemented(
          "StructuredColumns::createField currently relies"
          " on ATLAS_HAVE_TRANS for parallel fields", Here());
    }
    field::Field* field = field::Field::create<double>(
                    name, array::make_shape(grid_.npts(), levels));

    field->set_functionspace(*this);
    set_field_metadata(options,*field);
    return field;
#endif
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::gather(
    const field::FieldSet& local_fieldset,
    field::FieldSet& global_fieldset ) const
{
#ifdef ATLAS_HAVE_TRANS
    ASSERT(local_fieldset.size() == global_fieldset.size());

    for( size_t f=0; f<local_fieldset.size(); ++f )
    {
        const field::Field& loc = local_fieldset[f];
        field::Field& glb = global_fieldset[f];
        if( loc.datatype() != array::DataType::str<double>() )
        {
            std::stringstream err;
            err << "Cannot gather Structured field " << loc.name()
                << " of datatype " << loc.datatype().str() << ".";
            err << "Only " << array::DataType::str<double>() << " supported.";
            throw eckit::BadValue(err.str());
        }

        size_t root(0);
        glb.metadata().get("owner",root);
        std::vector<int> nto(1,root+1);
        if( loc.rank() > 1 )
        {
            nto.resize(loc.stride(0));
            for( size_t i=0; i<nto.size(); ++i )
            {
                nto[i] = root+1;
            }
        }
        trans_->gathgrid(nto.size(), nto.data(),
                         loc.data<double>(), glb.data<double>());
    }

#else
    eckit::NotImplemented("StructuredColumns::gather currently relies "
                          "on ATLAS_HAVE_TRANS", Here());
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void StructuredColumns::gather(
    const field::Field& local,
    field::Field& global) const
{
    field::FieldSet local_fields;
    field::FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields,global_fields);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::scatter(
    const field::FieldSet& global_fieldset,
    field::FieldSet& local_fieldset) const
{
#ifdef ATLAS_HAVE_TRANS
    ASSERT(local_fieldset.size() == global_fieldset.size());

    for( size_t f=0; f<local_fieldset.size(); ++f )
    {
        const field::Field& glb = global_fieldset[f];
        field::Field& loc = local_fieldset[f];
        if( loc.datatype() != array::DataType::str<double>() )
        {
            std::stringstream err;
            err << "Cannot scatter Structured field " << glb.name()
                << " of datatype " << glb.datatype().str() << ".";
            err << "Only " << array::DataType::str<double>() << " supported.";
            throw eckit::BadValue(err.str());
        }

        size_t root(0);
        glb.metadata().get("owner",root);
        std::vector<int> nfrom(1,root+1);
        if( loc.rank() > 1 )
        {
            nfrom.resize(loc.stride(0));
            for( size_t i=0; i<nfrom.size(); ++i)
            {
                nfrom[i] = root+1;
            }
        }
        trans_->distgrid(nfrom.size(), nfrom.data(),
                         glb.data<double>(), loc.data<double>());
        glb.metadata().broadcast(loc.metadata(),root);
        loc.metadata().set("global",false);
  }
#else
    eckit::NotImplemented("StructuredColumns::scatter currently relies "
                          "on ATLAS_HAVE_TRANS", Here());
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void StructuredColumns::scatter(
    const field::Field& global,
    field::Field& local) const
{
    field::FieldSet global_fields;
    field::FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}
// ----------------------------------------------------------------------------






// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double StructuredColumns::lat(
    size_t jlat) const
{
  return grid_.y(jlat+first_lat_);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double StructuredColumns::lon(
    size_t jlat,
    size_t jlon) const
{
  return grid_.x(jlon+first_lon_[jlat],jlat+first_lat_);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Checksum FieldSet
// ----------------------------------------------------------------------------
std::string StructuredColumns::checksum(
    const field::FieldSet& fieldset) const
{
    eckit::MD5 md5;
    NOTIMP;
    return md5;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Checksum Field
// ----------------------------------------------------------------------------
std::string StructuredColumns::checksum(
    const field::Field& field) const
{
    // FieldSet fieldset;
    // fieldset.add(field);
    // return checksum(fieldset);
    eckit::Log::warning() << "Only local checksum implemented" << std::endl;
    std::stringstream resultss;
    resultss << internals::checksum(field.data<double>(),field.size());
    return resultss.str();
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------
extern "C"
{

StructuredColumns* atlas__functionspace__StructuredColumns__new__grid (const grid::Grid* grid)
{
  ATLAS_ERROR_HANDLING(
    return new StructuredColumns(*grid);
  );
  return 0;
}

void atlas__functionspace__StructuredColumns__delete (StructuredColumns* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

field::Field* atlas__fs__StructuredColumns__create_field_name_kind (const StructuredColumns* This, const char* name, int kind, const eckit::Parametrisation* options)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),array::DataType(kind),*options);
  );
  return 0;
}

field::Field* atlas__fs__StructuredColumns__create_field_name_kind_lev (const StructuredColumns* This, const char* name, int kind, int levels, const eckit::Parametrisation* options)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),array::DataType(kind),levels,*options);
  );
  return 0;
}

void atlas__functionspace__StructuredColumns__gather (const StructuredColumns* This, const field::Field* local, field::Field* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->gather(*local,*global);
  );
}

void atlas__functionspace__StructuredColumns__scatter (const StructuredColumns* This, const field::Field* global, field::Field* local)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->scatter(*global,*local);
  );
}

void atlas__fs__StructuredColumns__checksum_fieldset(const StructuredColumns* This, const field::FieldSet* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__fs__StructuredColumns__checksum_field(const StructuredColumns* This, const field::Field* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

}
// ----------------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

