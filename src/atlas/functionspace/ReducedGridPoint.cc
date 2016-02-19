/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/utils/MD5.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/ReducedGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/ReducedGridPoint.h"
#include "atlas/private/Checksum.h"
#include "atlas/util/runtime/ErrorHandling.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/numerics/trans/Trans.h"
#endif

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
ReducedGridPoint::ReducedGridPoint(const Grid& grid) :
  FunctionSpace()
{
    grid_ = dynamic_cast<const grids::ReducedGrid*>(&grid);
    if (grid_ == NULL)
    {
      throw eckit::BadCast("Grid is not a grids::ReducedGrid type. "
                           "Cannot partition using IFS trans", Here());
    }

#ifdef ATLAS_HAVE_TRANS
    trans_ = new trans::Trans(*grid_);

    npts_ = trans_->ngptot();

    // Maximum number of global points in the longitude direction
    int nlonmax = dynamic_cast<const grids::ReducedGrid*>(&grid)->nlonmax();

    // Number of latitude bands
    int n_regions_NS = trans_->n_regions_NS();

    // Number of partitions per latitude band
    ArrayView<int,1> n_regions = trans_->n_regions();

    // First latitude of latitude band
    ArrayView<int,1> nfrstlat = trans_->nfrstlat();

    // First latitude of latitude band
    ArrayView<int,1> nlstlat = trans_->nlstlat();

    // Index of latitude partition (note that if a partition
    // has two regions on a latitude - the index increases
    // by one (2 numbers)
    ArrayView<int,1> nptrfrstlat = trans_->nptrfrstlat();

    // Starting longitudinal point per given latitude (ja)
    // Note that it is associated to nptrfrstlat
    ArrayView<int,2> nsta = trans_->nsta();

    // Number of longitudinal points per given latitude (ja)
    // Note that it is associated to nptrfrstlat
    ArrayView<int,2> nonl = trans_->nonl();

    // Total number of longitudes per latitude
    ArrayView<int,1> nloen = trans_->nloen();

    size_t proc(0);
    // Loop over number of latitude bands (ja)
    for (size_t ja = 0; ja < n_regions_NS; ++ja)
    {
        // Loop over number of longitude bands (jb)
        for (size_t jb = 0; jb < n_regions[ja]; ++jb)
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
    exit_outer_loop:


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
#else
    npts_ = grid_->npts();
    nlat_ = grid_->nlat();
    nlon_.resize(nlat_);
    for( size_t jlat=0; jlat<nlat_; ++jlat )
      nlon_[jlat] = grid_->nlon(jlat);
    first_lat_ = 0;
    first_lon_.resize(nlat_,0);
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------
ReducedGridPoint::~ReducedGridPoint()
{
#ifdef ATLAS_HAVE_TRANS
    delete trans_;
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
template <>
Field* ReducedGridPoint::createField<double>(
    const std::string& name) const
{
#ifdef ATLAS_HAVE_TRANS
    Field* field = Field::create<double>(name, make_shape(trans_->ngptot()));
    field->set_functionspace(this);
    return field;
#else
    eckit::NotImplemented("ReducedGridPoint::createField currently relies"
                          " on ATLAS_HAVE_TRANS", Here());

    Field* field = Field::create<double>(name, make_shape(grid_->npts()) );
    field->set_functionspace(this);
    return field;
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create Field with vertical levels
// ----------------------------------------------------------------------------
template <>
Field* ReducedGridPoint::createField<double>(
    const std::string& name,
    size_t levels) const
{
#ifdef ATLAS_HAVE_TRANS
    Field* field = Field::create<double>(
                    name, make_shape(trans_->ngptot(), levels));

    field->set_functionspace(this);
    field->set_levels(levels);
    return field;
#else
    eckit::NotImplemented("ReducedGridPoint::createField currently relies "
                          "on ATLAS_HAVE_TRANS", Here());

    Field* field = Field::create<double>(
                    name, make_shape(grid_->npts(), levels));

    field->set_functionspace(this);
    return field;
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create global Field
// ----------------------------------------------------------------------------
template <>
Field* ReducedGridPoint::createGlobalField<double>(
    const std::string& name) const
{
    Field* field = Field::create<double>(name, make_shape(grid_->npts()));
    field->set_functionspace(this);
    return field;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create global Field with vertical levels
// ----------------------------------------------------------------------------
template <>
Field* ReducedGridPoint::createGlobalField<double>(
    const std::string& name,
    size_t levels) const
{
    Field* field = Field::create<double>(
                    name, make_shape(grid_->npts(), levels));

    field->set_functionspace(this);
    return field;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void ReducedGridPoint::gather(
    const FieldSet& local_fieldset,
    FieldSet& global_fieldset ) const
{
#ifdef ATLAS_HAVE_TRANS
    ASSERT(local_fieldset.size() == global_fieldset.size());

    for( size_t f=0; f<local_fieldset.size(); ++f )
    {
        const Field& loc = local_fieldset[f];
        Field& glb = global_fieldset[f];
        if( loc.datatype() != DataType::str<double>() )
        {
            std::stringstream err;
            err << "Cannot gather ReducedGrid field " << loc.name()
                << " of datatype " << loc.datatype().str() << ".";
            err << "Only " << DataType::str<double>() << " supported.";
            throw eckit::BadValue(err.str());
        }

        std::vector<int> nto(1,1);
        if( loc.rank() > 1 )
        {
            nto.resize(loc.stride(0));
            for( size_t i=0; i<nto.size(); ++i )
            {
                nto[i] = 1;
            }
        }
        trans_->gathgrid(nto.size(), nto.data(),
                         loc.data<double>(), glb.data<double>());
    }

#else
    eckit::NotImplemented("ReducedGridPoint::gather currently relies "
                          "on ATLAS_HAVE_TRANS", Here());
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void ReducedGridPoint::gather(
    const Field& local,
    Field& global) const
{
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields,global_fields);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void ReducedGridPoint::scatter(
    const FieldSet& global_fieldset,
    FieldSet& local_fieldset) const
{
#ifdef ATLAS_HAVE_TRANS
    ASSERT(local_fieldset.size() == global_fieldset.size());

    for( size_t f=0; f<local_fieldset.size(); ++f )
    {
        const Field& glb = global_fieldset[f];
        Field& loc = local_fieldset[f];
        if( loc.datatype() != DataType::str<double>() )
        {
            std::stringstream err;
            err << "Cannot scatter ReducedGrid field " << glb.name()
                << " of datatype " << glb.datatype().str() << ".";
            err << "Only " << DataType::str<double>() << " supported.";
            throw eckit::BadValue(err.str());
        }

        std::vector<int> nfrom(1,1);
        if( loc.rank() > 1 )
        {
            nfrom.resize(loc.stride(0));
            for( size_t i=0; i<nfrom.size(); ++i)
            {
                nfrom[i] = 1;
            }
        }
        trans_->distgrid(nfrom.size(), nfrom.data(),
                         glb.data<double>(), loc.data<double>());
  }
#else
    eckit::NotImplemented("ReducedGridPoint::scatter currently relies "
                          "on ATLAS_HAVE_TRANS", Here());
#endif
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void ReducedGridPoint::scatter(
    const Field& global,
    Field& local) const
{
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}
// ----------------------------------------------------------------------------






// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double ReducedGridPoint::lat(
    size_t jlat) const
{
  return grid_->lat(jlat+first_lat_);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double ReducedGridPoint::lon(
    size_t jlat,
    size_t jlon) const
{
  return grid_->lon(jlat+first_lat_, jlon+first_lon_[jlat]);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Checksum FieldSet
// ----------------------------------------------------------------------------
std::string ReducedGridPoint::checksum(
    const FieldSet& fieldset) const
{
    eckit::MD5 md5;
    NOTIMP;
    return md5;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Checksum Field
// ----------------------------------------------------------------------------
std::string ReducedGridPoint::checksum(
    const Field& field) const
{
    // FieldSet fieldset;
    // fieldset.add(field);
    // return checksum(fieldset);
    eckit::Log::warning() << "Only local checksum implemented" << std::endl;
    std::stringstream resultss;
    resultss << atlas::checksum(field.data<double>(),field.size());
    return resultss.str();
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------
extern "C"
{

ReducedGridPoint* atlas__functionspace__ReducedGridPoint__new__grid (const Grid* grid)
{
  ATLAS_ERROR_HANDLING(
    return new ReducedGridPoint(*grid);
  );
  return 0;
}

void atlas__functionspace__ReducedGridPoint__delete (ReducedGridPoint* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Field* atlas__functionspace__ReducedGridPoint__create_field (const ReducedGridPoint* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoint__create_field_lev (const ReducedGridPoint* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name),levels);
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoint__create_gfield (const ReducedGridPoint* This, const char* name)
{
  ATLAS_ERROR_HANDLING (
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoint__create_gfield_lev (const ReducedGridPoint* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name),levels);
  );
  return 0;
}

void atlas__functionspace__ReducedGridPoint__gather (const ReducedGridPoint* This, const Field* local, Field* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->gather(*local,*global);
  );
}

void atlas__functionspace__ReducedGridPoint__scatter (const ReducedGridPoint* This, const Field* global, Field* local)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->scatter(*global,*local);
  );
}

void atlas__fs__ReducedGridPoint__checksum_fieldset(const ReducedGridPoint* This, const FieldSet* fieldset, char* &checksum, int &size, int &allocated)
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

void atlas__fs__ReducedGridPoint__checksum_field(const ReducedGridPoint* This, const Field* field, char* &checksum, int &size, int &allocated)
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

