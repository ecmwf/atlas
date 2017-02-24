/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/Structured.h"

#include <algorithm>
#include <limits>
#include "atlas/runtime/ErrorHandling.h"
#include "eckit/geometry/Point3.h"
#include "atlas/internals/Debug.h"

using eckit::geometry::Point3;
using eckit::geometry::lonlat_to_3d;

namespace atlas {
namespace grid {


Structured* Structured::create(const util::Config& p) {
    Structured* grid = dynamic_cast<Structured*>(Grid::create(p));
    if (!grid)
        throw eckit::BadParameter("Grid is not a reduced grid", Here());
    return grid;

}


Structured* Structured::create(const std::string& uid) {
    Structured* grid = dynamic_cast<Structured*>( Grid::create(uid) );
    if (!grid)
        throw eckit::BadParameter("Grid "+uid+" is not a reduced grid",Here());
    return grid;
}


std::string Structured::className() {
    return "atlas.grid.Structured";
}


std::string Structured::grid_type_str() {
    return "structured";
}

std::string Structured::shortName() const {
  return "structured";

}

Structured::Structured() :
    Grid(),
    N_(0) {
}


Structured::~Structured() {
}

/*
void Structured::setup(
    const size_t ny,
    const double y[],
    const long nx[],
    const double xmin[],
    const double xmax[] ) {
    ASSERT(ny>1);  // can't have a grid with just one latitude

    nx_  .assign(nx,   nx   + ny);
    y_   .assign(y,    y    + ny);
    xmin_.assign(xmin, xmin + ny);
    xmax_.assign(xmax, xmax + ny);
    npts_ = static_cast<size_t>(std::accumulate(nx_.begin(), nx_.end(), 0));

    dx_.resize(ny);
    nxmin_ = nxmax_ = static_cast<size_t>(nx_[0]);

    for (size_t j = 0; j < ny; ++j) {
        dx_[j] = (xmax_[j]-xmin_[j])/double(nx_[j]-1);
        nxmin_ = std::min(static_cast<size_t>(nx_[j]),nxmin_);
        nxmax_ = std::max(static_cast<size_t>(nx_[j]),nxmax_);
    }

    compute_true_periodicity();
}



void Structured::setup_cropped(const size_t ny, const double y[], const long nx[], const double xmin[], const double xmax[], const domain::Domain& dom ) {
    ASSERT(ny>0);

    std::vector<double> dom_y;    dom_y.  reserve(ny);
    std::vector<long>   dom_nx;   dom_nx.    reserve(ny);
    std::vector<double> dom_xmin; dom_xmin.reserve(ny);
    std::vector<double> dom_xmax; dom_xmax.reserve(ny);
    const double tol = 1.e-6;
    size_t dom_ny = 0;
    const double d_xmin = dom.xmin();
    const double d_xmax = dom.xmax();
    const double d_ymin = dom.ymin();
    const double d_ymax = dom.ymax();
    double dx;
    for( size_t j=0; j<ny; ++j )
    {
        if( y[j]-tol < d_ymax && y[j]+tol > d_ymin )
        {
            ++dom_ny;
            const double _y = y[j];
            double _xmin = xmin[j];
            double _xmax = xmax[j];
            if( isPeriodicX() ) {  // periodic:      nx = number of divisions
              dx = (d_xmax-d_xmin)/double(nx[j]);
            } else {            // not periodic:  nx = number of points
              if( _xmin < _xmax ) {
                dx = (_xmax-_xmin)/double(nx[j]-1l);
              } else {
                dx = (d_xmax-d_xmin)/double(nx[j]-1l);
              }
            }


// Flaw: Assumed that original grid has increments that start from xmin=0.0
//       xmin of original grid does not even have to be multiple of dx

            long _nx(0);
            for( long i=0; i<nx[j]; ++i )
            {
                const double x = dx*i;
                if( x+tol > d_xmin && x-tol < d_xmax )
                {
                    _xmin = std::min(_xmin,x);
                    _nx++;
                }
            }
            _xmax = _xmin+_nx*dx;
            dom_y     .push_back(y[j]);
            dom_nx    .push_back(_nx);
            dom_xmin.push_back(_xmin);
            dom_xmax.push_back(_xmax);
        }
    }
    setup(dom_ny,dom_y.data(),dom_nx.data(),dom_xmin.data(),dom_xmax.data());
}
*/


void Structured::setup( 
    spacing::Spacing*          yspace,
    const std::vector<long>&   nx,
    const std::vector<double>& xmin,
    const std::vector<double>& xmax,
    const std::vector<double>& dx )
{
  yspace_ = std::unique_ptr<spacing::Spacing>(yspace);
  size_t ny = yspace_->size();

  ASSERT(nx.size() == ny);
  ASSERT(dx.size() == ny);
  ASSERT(xmin.size() == ny);
  ASSERT(xmax.size() == ny);

  y_.assign(yspace->begin(),yspace->end());
  nx_.assign(nx.begin(),nx.end());
  dx_.assign(dx.begin(),dx.end());
  xmin_.assign(xmin.begin(),xmin.end());
  xmax_.assign(xmax.begin(),xmax.end());
  
  nxmin_ = nxmax_ = static_cast<size_t>(nx_.front());

  for (size_t j=1; j<ny; ++j) {
      nxmin_ = std::min(static_cast<size_t>(nx_[j]),nxmin_);
      nxmax_ = std::max(static_cast<size_t>(nx_[j]),nxmax_);
  }

  npts_ = static_cast<size_t>(std::accumulate(nx_.begin(), nx_.end(), 0));

  compute_true_periodicity();
}

void Structured::setup( 
    spacing::Spacing* yspace,
    const long&       nx,
    const double&     xmin,
    const double&     xmax,
    const double&     dx )
{
  yspace_ = std::unique_ptr<spacing::Spacing>(yspace);
  size_t ny = yspace_->size();

  y_.assign(yspace->begin(),yspace->end());
  nx_.assign(ny, nx);
  dx_.assign(ny, dx);
  xmin_.assign(ny, xmin);
  xmax_.assign(ny, xmax);
  
  nxmin_ = nxmax_ = static_cast<size_t>(nx_.front());

  npts_ = static_cast<size_t>(std::accumulate(nx_.begin(), nx_.end(), 0));

  compute_true_periodicity();
}

namespace {
  struct EarthCentred {
    EarthCentred(Structured& _grid) : grid(_grid) {}
    Point3 operator()(double ll[]) {
      grid.projection().xy2lonlat(ll);
      double p[3];
      double r=1., h=0.;
      lonlat_to_3d(ll,p,r,h);
      return Point3( p[0], p[1], p[2] );
    }
    Point3 operator()(std::array<double,2> ll) { return operator()(ll.data()); }
    Structured& grid;
  };
}

void Structured::compute_true_periodicity() {

  if( projection_->isStrictlyRegional() ) {
    periodic_x_ = false;

  } else if( domain_->isGlobal() ) {
    periodic_x_ = true;

  } else {
    // domain could be zonal band

    size_t j = ny()/2;
    if( (nx_[j]-1) * dx_[j] == xmax_[j] ) {
      periodic_x_ = false; // This would lead to duplicated points
    } else {
      // High chance to be periodic. Check anyway.
      EarthCentred compute_earth_centred(*this);

      Point3 Pxmin = compute_earth_centred( { xmin_[j], y_[j] } );
      Point3 Pxmax = compute_earth_centred( { xmax_[j], y_[j] } );

      periodic_x_  = points_equal( Pxmin, Pxmax );      
    }

  }

}

void Structured::lonlat( std::vector<Point>& pts ) const {
    pts.resize(npts());

    for(size_t jlat=0, c=0; jlat<nlat(); ++jlat) {
        const double y = lat(jlat);
        for(size_t jlon=0; jlon<nlon(jlat); ++jlon) {
            pts[c++].assign(lon(jlat,jlon),y);
        }
    }
}


size_t Structured::copyLonLatMemory(double* pts, size_t size) const {
    size_t sizePts = 2*npts();
    ASSERT(size >= sizePts);

    for(size_t jlat=0, c=0; jlat<nlat(); ++jlat ) {
        const double y = lat(jlat);
        for( size_t jlon=0; jlon<nlon(jlat); ++jlon ) {
            pts[c++] = lon(jlat,jlon);
            pts[c++] = y;
        }
    }
    return sizePts;
}


void Structured::print(std::ostream& os) const {
    os << "Structured(Name:" << shortName() << ")";
}


void Structured::hash(eckit::MD5& md5) const {
    // Through inheritance the grid_type_str() might differ while still being same grid
    //md5.add(grid_type_str());

    md5.add(latitudes().data(), sizeof(double)*latitudes().size());
    md5.add(pl().data(), sizeof(long)*nlat());

    // also add lonmin and lonmax
    md5.add(xmin_.data(), sizeof(double)*xmin_.size());
    md5.add(xmax_.data(), sizeof(double)*xmax_.size());

    // also add projection information
    eckit::Properties prop;
    std::ostringstream s;
    s << projection().spec();
    prop.set("projection",s.str());
    prop.hash(md5);
}

eckit::Properties Structured::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specific specs
    grid_spec.set("nlat",nlat());
    grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));
    grid_spec.set("pl",eckit::makeVectorValue(pl()));
    grid_spec.set("lonmin",eckit::makeVectorValue(xmin_));
    grid_spec.set("lonmax",eckit::makeVectorValue(xmax_));

    return grid_spec;
}


// --------------------------------------------------------------------


extern "C" {


    size_t atlas__grid__Structured__N(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->N();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlat(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlat();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlon(Structured* This, size_t jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlon(jlat);
        );
        return 0;
    }


    void atlas__grid__Structured__pl(Structured* This, const long* &nlons, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            nlons = This->pl().data();
            size  = This->pl().size();
        );
    }


    size_t atlas__grid__Structured__nlonmax(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlonmax();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlonmin(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlonmin();
        );
        return 0;
    }


    size_t atlas__grid__Structured__npts(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->npts();
        );
        return 0;
    }


    double atlas__grid__Structured__lat(Structured* This,size_t jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->lat(jlat);
        );
        return 0.;
    }


    double atlas__grid__Structured__lon(Structured* This,size_t jlat,size_t jlon) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->lon(jlat, jlon);
        );
        return 0.;
    }


    void atlas__grid__Structured__lonlat(Structured* This, size_t jlat, size_t jlon, double crd[]) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            This->lonlat(jlat, jlon, crd);
        );
    }


    void atlas__grid__Structured__latitudes(Structured* This, const double* &lat, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            lat  = This->latitudes().data();
            size = This->latitudes().size();
        );
    }


    int atlas__grid__Structured__reduced(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->reduced();
        );
        return 1;
    }


    Structured* atlas__grid__Structured(char* identifier) {
        ATLAS_ERROR_HANDLING(
            ASSERT( identifier );
            return Structured::create( std::string(identifier) );
        );
        return 0;
    }


    Structured* atlas__grid__Structured__config(util::Config* conf) {
        ATLAS_ERROR_HANDLING(
            ASSERT( conf );
            return Structured::create(*conf);
        );
        return 0;
    }


    void atlas__grid__Structured__delete(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
        );
        delete This;
    }


}


}  // namespace grid
}  // namespace atlas

