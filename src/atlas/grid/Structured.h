/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef Structured_h
#define Structured_h

#include "eckit/memory/Builder.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/Config.h"

#include "atlas/grid/spacing/Spacing.h"


namespace atlas {
namespace grid {


/**
 * @brief Structured Grid
 *
 * This class is a base class for all grids that can be described by
 * constant latitudes with a uniform distribution of points per latitude
 * in zonal direction.
 * This means any full grid and reduced grid, both regular, gaussian or other
 * such distribution can be represented with this class
 */
class Structured : public Grid {

public:

    typedef eckit::SharedPtr<Structured> Ptr;

    static Structured* create(const util::Config&);

    static Structured* create(const std::string& shortname);

public:

    static std::string className();

    static std::string grid_type_str();

public:

    Structured();

    virtual ~Structured();

    virtual size_t npts() const {
        return npts_;
    }

    virtual void lonlat(std::vector<Point>&) const;

    //virtual std::string gridType() const {
    //    return grid_type_;
     //}

    virtual eckit::Properties spec() const;

    /**
     * Human readable name
     * @note: may not be unique, such as when reduced Gaussian grids have the same N numbers but different distribution of latitude points
     */
    virtual std::string shortName() const;
    virtual std::string gridType() const { return "structured"; }


    virtual std::string getOptimalMeshGenerator() const {
        return "Structured";
    }

    virtual size_t N() const {
        return N_;
    }

    inline size_t nlat() const {
        return y_.size();
    }

    inline size_t nlon( size_t jlat ) const {
        return static_cast<size_t>(nx_[jlat]);
    }

    inline size_t nlonmax() const {
        return nxmax_;
    }

    inline size_t nlonmin() const {
        return nxmin_;
    }

    inline const std::vector<long>& pl() const {
        return nx_;
    }

    inline const std::vector<double>& latitudes() const {
        return y_;
    }

    inline double lon( const size_t jlat, const size_t jlon ) const {
        return xmin_[jlat] + static_cast<double>(jlon) * dx_[jlat];
    }

    inline double lat( const size_t jlat ) const {
        return y_[jlat];
    }

    inline void lonlat( const size_t jlat, const size_t jlon, double crd[] ) const {
        crd[0] = lon(jlat,jlon);
        crd[1] = lat(jlat);
    }

    void geoLonlat(const size_t jlon, const size_t jlat, eckit::geometry::LLPoint2 &Pll) const {
      // in grid coordinates
      double crd[2];
      lonlat(jlat,jlon,crd);
      // convert to geographic coordinates
      projection_->coords2lonlat(crd);
      Pll.assign(crd[0],crd[1]);
    }

    inline bool reduced() const {
        return nlonmax() != nlonmin();
    }

protected: // methods

    virtual size_t copyLonLatMemory(double* pts, size_t size) const;

    virtual void print(std::ostream&) const;

    /// Hash of the PL array
    virtual void hash(eckit::MD5&) const;

    void setup(const size_t nlat, const double lats[], const long pl[], const double lonmin[], const double lonmax[]);

    void setup_cropped(const size_t ny, const double y[], const long nx[], const double xmin[], const double xmax[], const domain::Domain& dom);

protected:

    /// Number of latitudes in hemisphere
    size_t N_;

    // TODO: document
    size_t nxmin_;

    // TODO: document
    size_t nxmax_;

    /// Total number of unique points in the grid
    size_t npts_;

    /// Latitude values
    std::vector<double> y_;

    // TODO: remove, only to instantiate leaf classes
    //std::string grid_type_;

    // TODO: remove, only to instantiate leaf classes
    //std::string shortName_;

    /// Number of points per latitude
    std::vector<long> nx_;

    /// Value of minimum longitude per latitude [default=0]
    std::vector<double> xmin_;

    /// Value of maximum longitude per latitude [default=0]
    std::vector<double> xmax_;

    /// Value of longitude increment
    std::vector<double> dx_;

};


#define CONFIG util::Config
extern "C"
{
    void atlas__grid__Structured__delete(Structured* This);
    Structured* atlas__grid__Structured(char* identifier);
    Structured* atlas__grid__Structured__config(CONFIG* conf);
    Structured* atlas__grid__CustomStructured_int(size_t nlat, double lat[], int nlon[]);
    Structured* atlas__grid__CustomStructured_long(size_t nlat, double lat[], long nlon[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(size_t nlat, double lat[], int nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(size_t nlat, double lat[], long nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__regular__RegularGaussian(size_t N);
    Structured* atlas__grid__reduced__ReducedGaussian_int(size_t N, int nlon[]);
    Structured* atlas__grid__reduced__ReducedGaussian_long(size_t N, long nlon[]);
    Structured* atlas__grid__regular__RegularLonLat(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLonLat(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLon(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLat(size_t nlon, size_t nlat);

    void   atlas__grid__Structured__pl        (Structured* This, const long* &pl, size_t &size);
    size_t atlas__grid__Structured__N         (Structured* This);
    size_t atlas__grid__Structured__nlat      (Structured* This);
    size_t atlas__grid__Structured__nlon      (Structured* This, size_t jlat);
    size_t atlas__grid__Structured__nlonmin   (Structured* This);
    size_t atlas__grid__Structured__nlonmax   (Structured* This);
    size_t atlas__grid__Structured__npts      (Structured* This);
    double atlas__grid__Structured__lat       (Structured* This, size_t jlat);
    double atlas__grid__Structured__lon       (Structured* This, size_t jlat, size_t jlon);
    void   atlas__grid__Structured__lonlat    (Structured* This, size_t jlat, size_t jlon, double crd[]);
    void   atlas__grid__Structured__latitudes (Structured* This, const double* &lats, size_t &size);
    int    atlas__grid__Structured__reduced   (Structured* This);
}
#undef CONFIG


}  // namespace grid
}  // namespace atlas


#endif
