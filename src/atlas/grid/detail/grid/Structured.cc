/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "Structured.h"

#include <algorithm>
#include <limits>
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/Point.h"
#include "atlas/internals/Debug.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/grid/detail/grid/GridBuilder.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {


std::string Structured::static_type() {
    return "structured";
}

std::string Structured::name() const {
  return name_;
}

Structured::Structured( Projection p, XSpace* xspace, YSpace yspace, Domain domain ):
    Structured( Structured::static_type(), p, xspace, yspace, domain ) {
}

Structured::Structured( const std::string& name, Projection projection, XSpace* xspace, YSpace yspace, Domain domain):
  Grid(),
  name_(name) {
  // Copry members
  if( projection )
    projection_ = projection;
  else
    projection_ = Projection();
  xspace_ = std::unique_ptr<XSpace>( xspace );
  yspace_ = yspace;
  domain_ = domain;
  ASSERT( xspace_->ny == yspace_.size() );

  y_.assign(yspace_.begin(),yspace_.end());
  nx_    = xspace_->nx;
  dx_    = xspace_->dx;
  xmin_  = xspace_->xmin;
  xmax_  = xspace_->xmax;

  // Further setup
  nxmin_ = nxmax_ = static_cast<size_t>(nx_.front());
  size_t ny = xspace_->ny;
  for (size_t j=1; j<ny; ++j) {
      nxmin_ = std::min(static_cast<size_t>(nx_[j]),nxmin_);
      nxmax_ = std::max(static_cast<size_t>(nx_[j]),nxmax_);
  }
  npts_ = size_t(std::accumulate(nx_.begin(), nx_.end(), 0));
  compute_true_periodicity();

  Log::info() << Here() << std::endl;;
  if( periodic() ) {
    Log::info() << Here() << std::endl;;
    if( yspace.max() - yspace.min() == 180. ) {
      Log::info() << Here() << std::endl;;
      domain_ = Domain(Grid::Config("type","global"));
    }
    else {
      Log::info() << Here() << std::endl;;
      Grid::Config dom;
      dom.set("type","zonal_band");
      dom.set("ymin",yspace.min());
      dom.set("ymax",yspace.max());
      domain_ = Domain(dom);
    }
  }
}


Structured::~Structured() {
}

Structured::XSpace::XSpace(long _ny) :
  ny(_ny) {
  nx.  reserve(ny);
  xmin.reserve(ny);
  xmax.reserve(ny);
  dx  .reserve(ny);
}

Structured::XSpace::XSpace( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint ) :
  ny(N.size()),
  nx(N),
  xmin(ny,interval[0]),
  xmax(ny,interval[1]),
  dx(ny) {
  nxmin = std::numeric_limits<size_t>::max();
  nxmax = 0;
  double length = interval[1] - interval[0];
  for( size_t j=0; j<ny; ++j ) {
    nxmin = std::min( nxmin, size_t(nx[j]) );
    nxmax = std::max( nxmax, size_t(nx[j]) );
    dx[j] = endpoint ? length/double(nx[j]-1) : length/double(nx[j]);
  }
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
    const YSpace&              yspace,
    const std::vector<long>&   nx,
    const std::vector<double>& xmin,
    const std::vector<double>& xmax,
    const std::vector<double>& dx )
{
  yspace_ = yspace;
  size_t ny = yspace_.size();

  ASSERT(nx.size() == ny);
  ASSERT(dx.size() == ny);
  ASSERT(xmin.size() == ny);
  ASSERT(xmax.size() == ny);

  y_.assign(yspace_.begin(),yspace_.end());
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

/*
void Structured::setup(
    const YSpace&     yspace,
    const long&       nx,
    const double&     xmin,
    const double&     xmax,
    const double&     dx )
{
  yspace_ = yspace;
  size_t ny = yspace_.size();

  y_.assign(yspace_.begin(),yspace_.end());
  nx_.assign(ny, nx);
  dx_.assign(ny, dx);
  xmin_.assign(ny, xmin);
  xmax_.assign(ny, xmax);

  nxmin_ = nxmax_ = static_cast<size_t>(nx_.front());

  npts_ = static_cast<size_t>(std::accumulate(nx_.begin(), nx_.end(), 0));

  compute_true_periodicity();
}
*/

namespace {
  struct EarthCentred {
    EarthCentred(Structured& _grid) : grid(_grid) {}
    PointXYZ operator()(const PointXY& xy) {
      PointLonLat lonlat = grid.projection().lonlat(xy);
      return lonlat_to_geocentric(lonlat,radius);
    }
    Structured& grid;
    double radius={1.};
  };
}

void Structured::compute_true_periodicity() {

  if( projection_.isStrictlyRegional() ) {
    periodic_x_ = false;

  } else if( domain_.global() ) {
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


void Structured::print(std::ostream& os) const {
    os << "Structured(Name:" << name() << ")";
}

std::string Structured::type() const {
  if( not type_.empty() ) {
    type_ += std::string(reduced()?"reduced":"regular");
  }
  return type_;
}

void Structured::hash(eckit::MD5& md5) const {
    // Through inheritance the static_type() might differ while still being same grid
    //md5.add(static_type());

    md5.add(y().data(), sizeof(double)*y().size());
    md5.add(nx().data(), sizeof(long)*ny());

    // also add lonmin and lonmax
    md5.add(xmin_.data(), sizeof(double)*xmin_.size());
    md5.add(xmax_.data(), sizeof(double)*xmax_.size());

    // also add projection information
    Grid::Spec prop;
    std::ostringstream s;
    s << projection().spec();
    prop.set("projection",s.str());
    prop.hash(md5);
}

Grid::Spec Structured::spec() const {
    Grid::Spec grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specific specs
    grid_spec.set("yspace",yspace().spec());
    // grid_spec.set("y",eckit::makeVectorValue(latitudes()));
    grid_spec.set("nx[]",eckit::makeVectorValue(nx()));
    grid_spec.set("xmin",eckit::makeVectorValue(xmin_));
    grid_spec.set("xmax",eckit::makeVectorValue(xmax_));

    return grid_spec;
}


// --------------------------------------------------------------------

namespace { // anonymous

static class structured : public GridBuilder {

public:

    structured(): GridBuilder( "structured" ){}

    virtual void print(std::ostream& os) const {
        os << std::left << std::setw(20) << " " << "Structured grid";
    }

    virtual const atlas::grid::Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
        throw eckit::NotImplemented( "Cannot create structured grid from name", Here() );
    }
  
    virtual const atlas::grid::Grid::grid_t* create( const Grid::Config& config ) const {

        Projection projection;
        Spacing    yspace;
        Domain     domain;

        Grid::Config config_proj;
        if( config.get("projection",config_proj) )
            projection = Projection(config_proj);

        Grid::Config config_yspace;
        if( not config.get("yspace",config_yspace) )
            throw eckit::BadParameter("yspace missing in configuration");
        yspace = Spacing(config_yspace);

        const size_t ny = yspace.size();

        StructuredGrid::grid_t::XSpace *X = new StructuredGrid::grid_t::XSpace(ny);

        double dom_xmin =  std::numeric_limits<double>::max();
        double dom_xmax = -std::numeric_limits<double>::max();

        std::vector<Grid::Config> config_xspace_list;
        if( config.get("xspace[]",config_xspace_list) ) {

            ASSERT( config_xspace_list.size() == ny );
            std::string xspace_type;

            for( size_t j=0; j<ny; ++j ) {
                config_xspace_list[j].get("type",xspace_type);
                ASSERT( xspace_type == "linear" );
                spacing::LinearSpacing::Params xspace( config_xspace_list[j] );
                X->xmin.push_back(xspace.start);
                X->xmax.push_back(xspace.end);
                X->nx.push_back(xspace.N);
                X->dx.push_back(xspace.step);
                dom_xmin = std::min(dom_xmin,xspace.start);
                dom_xmax = std::max(dom_xmax,xspace.end);
            }

        } else {

            Grid::Config config_xspace;
            if( not config.get("xspace",config_xspace) )
              throw eckit::BadParameter("xspace missing in configuration");

            std::string xspace_type;
            config_xspace.get("type",xspace_type);
            ASSERT( xspace_type == "linear" );

            std::vector<long>   v_N;
            std::vector<double> v_start;
            std::vector<double> v_end;
            std::vector<double> v_length;
            config_xspace.get("N[]",      v_N     );
            config_xspace.get("start[]",  v_start );
            config_xspace.get("end[]",    v_end   );
            config_xspace.get("length[]", v_length);

            if( not v_N.     empty() ) ASSERT(v_N.     size() == ny);
            if( not v_start. empty() ) ASSERT(v_start. size() == ny);
            if( not v_end.   empty() ) ASSERT(v_end.   size() == ny);
            if( not v_length.empty() ) ASSERT(v_length.size() == ny);

            for( size_t j=0; j<ny; ++j ) {
              if( not v_N.     empty() ) config_xspace.set("N",     v_N[j]);
              if( not v_start. empty() ) config_xspace.set("start", v_start[j]);
              if( not v_end.   empty() ) config_xspace.set("end",   v_end[j]);
              if( not v_length.empty() ) config_xspace.set("length",v_length[j]);
              spacing::LinearSpacing::Params xspace( config_xspace );
              X->xmin.push_back(xspace.start);
              X->xmax.push_back(xspace.end);
              X->nx.push_back(xspace.N);
              X->dx.push_back(xspace.step);
              dom_xmin = std::min(dom_xmin,xspace.start);
              dom_xmax = std::max(dom_xmax,xspace.end);
            }
        }

        Grid::Config config_domain;
        if( config.get("domain",config_domain) ) {
            domain = Domain(config_domain);
        } else {
            config_domain.set("type","rectangular");
            config_domain.set("ymin",yspace.min());
            config_domain.set("ymax",yspace.max());
            config_domain.set("xmin",dom_xmin);
            config_domain.set("xmax",dom_xmax);
            config_domain.set("units",projection.units());
            domain = Domain(config_domain);
        }

        return new StructuredGrid::grid_t(projection, X, yspace, domain );
    }

} structured_;

} // anonymous namespace


// --------------------------------------------------------------------


extern "C" {


    long atlas__grid__Structured__nlat(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->ny();
        );
        return 0;
    }


    long atlas__grid__Structured__nlon(Structured* This, long jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nx(jlat);
        );
        return 0;
    }


    void atlas__grid__Structured__pl(Structured* This, const long* &nlons, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            nlons = This->nx().data();
            size  = This->nx().size();
        );
    }


    long atlas__grid__Structured__nlonmax(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nxmax();
        );
        return 0;
    }


    long atlas__grid__Structured__nlonmin(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nxmin();
        );
        return 0;
    }


    long atlas__grid__Structured__npts(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->npts();
        );
        return 0;
    }


    double atlas__grid__Structured__lat(Structured* This, long jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->y(jlat);
        );
        return 0.;
    }


    double atlas__grid__Structured__lon( Structured* This, long jlat, long jlon ) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->x(jlon, jlat);
        );
        return 0.;
    }


    void atlas__grid__Structured__lonlat( Structured* This, long jlat, long jlon, double crd[] ) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            This->xy(jlon, jlat, crd);
        );
    }


    void atlas__grid__Structured__latitudes( Structured* This, const double* &lat, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            lat  = This->y().data();
            size = This->y().size();
        );
    }


    int atlas__grid__Structured__reduced(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->reduced();
        );
        return 1;
    }


    const Structured* atlas__grid__Structured(char* identifier) {
        ATLAS_ERROR_HANDLING(
            ASSERT( identifier );
            const Structured* grid = dynamic_cast<const Structured*>( Grid::create( std::string(identifier) ) );
            ASSERT( grid );
            return grid;
        );
        return 0;
    }


    const Structured* atlas__grid__Structured__config(util::Config* conf) {
        ATLAS_ERROR_HANDLING(
            ASSERT( conf );
            const Structured* grid = dynamic_cast<const Structured*>( Grid::create( *conf ) );
            ASSERT( grid );
            return grid;
        );
        return 0;
    }


    void atlas__grid__Structured__delete(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
        );
        delete This;
    }


    Structured* atlas__grid__regular__RegularGaussian(long N) {
      NOTIMP;
    }
    Structured* atlas__grid__regular__RegularLonLat(long nlon, long nlat) {
      NOTIMP;
    }
    Structured* atlas__grid__regular__ShiftedLonLat(long nlon, long nlat) {
      NOTIMP;
    }
    Structured* atlas__grid__regular__ShiftedLon(long nlon, long nlat) {
      NOTIMP;
    }
    Structured* atlas__grid__regular__ShiftedLat(long nlon, long nlat) {
      NOTIMP;
    }


    long atlas__grid__Gaussian__N(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            GaussianGrid gaussian(This);
            ASSERT( gaussian );
            return gaussian.N();
        );
        return 0;
    }
    

    Structured* atlas__grid__CustomStructured_int(long nlat, double lats[], int pl[]) {
      NOTIMP;
        // std::vector<long> pl_vector;
        // pl_vector.assign(pl,pl+nlat);
        // return new CustomStructured(nlat, lats, pl_vector.data());
    }


    Structured* atlas__grid__CustomStructured_long(long nlat, double lats[], long pl[]) {
      NOTIMP;
        // return new CustomStructured(nlat, lats, pl);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(long nlat, double lats[], int pl[], double lonmin[], double lonmax[]) {
      NOTIMP;
        // std::vector<long> pl_vector;
        // pl_vector.assign(pl, pl+nlat);
        // return new CustomStructured(nlat, lats, pl_vector.data(), lonmin, lonmax);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(long nlat, double lats[], long pl[], double lonmin[], double lonmax[]) {
        NOTIMP;
        //return new CustomStructured(nlat, lats, pl, lonmin, lonmax);
    }




}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

