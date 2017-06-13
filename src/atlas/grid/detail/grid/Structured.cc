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
#include "eckit/types/FloatCompare.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/Point.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"
#include "atlas/domain/Domain.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/runtime/Log.h"

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

Structured::Structured( XSpace xspace, YSpace yspace, Projection p, Domain domain ):
    Structured( Structured::static_type(), xspace, yspace, p, domain ) {
}

Structured::Structured( const std::string& name, XSpace xspace, YSpace yspace, Projection projection, Domain domain):
  Grid(),
  name_(name),
  xspace_(xspace),
  yspace_(yspace) {
  // Copry members
  if( projection )
    projection_ = projection;
  else
    projection_ = Projection();

  y_.assign(yspace_.begin(),yspace_.end());
  size_t ny = y_.size();

  if( xspace_.ny() == 1 && yspace_.size() > 1 ) {
    nx_   .resize( ny, xspace_.nx()[0]   );
    dx_   .resize( ny, xspace_.dx()[0]   );
    xmin_ .resize( ny, xspace_.xmin()[0] );
    xmax_ .resize( ny, xspace_.xmax()[0] );
  } else {
    nx_    = xspace_.nx();
    dx_    = xspace_.dx();
    xmin_  = xspace_.xmin();
    xmax_  = xspace_.xmax();
  }

  ASSERT( nx_.size() == ny );

  // Further setup
  nxmin_ = nxmax_ = static_cast<size_t>(nx_.front());
  for (size_t j=1; j<ny; ++j) {
      nxmin_ = std::min(static_cast<size_t>(nx_[j]),nxmin_);
      nxmax_ = std::max(static_cast<size_t>(nx_[j]),nxmax_);
  }
  npts_ = size_t(std::accumulate(nx_.begin(), nx_.end(), 0));


  if( not domain.empty() ) {
    crop( domain );
  }

  computeTruePeriodicity();

  if( domain.global() )
    domain_ = Domain( Grid::Config("type","global") );
  else
    computeDomain();
}

void Structured::computeDomain() {
  if( periodic() ) {
    if( yspace().max() - yspace().min() == 180. ) {
      domain_ = Domain( Grid::Config("type","global") );
    }
    else {
      Grid::Config config;
      config.set("type","zonal_band");
      config.set("ymin",yspace().min());
      config.set("ymax",yspace().max());
      domain_ = Domain(config);
    }
  } else if( domain_.empty() ) {
    Grid::Config config;
    config.set("type","rectangular");
    config.set("xmin",xmin_[0]);
    config.set("xmax",xmax_[0]);
    config.set("ymin",yspace().min());
    config.set("ymax",yspace().max());
    config.set("units",projection_.units());
    domain_ = Domain(config);
  }
}

Structured::~Structured() {
}

Structured::XSpace::XSpace() :
    impl_( nullptr ){
}

Structured::XSpace::XSpace( const XSpace& xspace ) :
    impl_( xspace.impl_ ){
}

Structured::XSpace::XSpace( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint ) :
    impl_( new Implementation(interval,N,endpoint) ) {
}

Structured::XSpace::XSpace( const Spacing& spacing ) :
    impl_( new Implementation(spacing) ) {
}

Structured::XSpace::XSpace( const Config& config ) :
    impl_( new Implementation(config) ) {
}

Structured::XSpace::XSpace( const std::vector<Config>& config ) :
    impl_( new Implementation(config) ) {
}

Structured::XSpace::Implementation::Implementation( const Config& config ) {

  Config config_xspace(config);

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

  long ny =  std::max( v_N.     size(),
             std::max( v_start. size(),
             std::max( v_end.   size(),
             std::max( v_length.size(),
             1ul ))));
  reserve(ny);

  if( not v_N.     empty() ) ASSERT(v_N.     size() == ny);
  if( not v_start. empty() ) ASSERT(v_start. size() == ny);
  if( not v_end.   empty() ) ASSERT(v_end.   size() == ny);
  if( not v_length.empty() ) ASSERT(v_length.size() == ny);

  nxmin_ = std::numeric_limits<size_t>::max();
  nxmax_ = 0;

  for( size_t j=0; j<ny; ++j ) {
    if( not v_N.     empty() ) config_xspace.set("N",     v_N[j]);
    if( not v_start. empty() ) config_xspace.set("start", v_start[j]);
    if( not v_end.   empty() ) config_xspace.set("end",   v_end[j]);
    if( not v_length.empty() ) config_xspace.set("length",v_length[j]);
    spacing::LinearSpacing::Params xspace( config_xspace );
    xmin_.push_back(xspace.start);
    xmax_.push_back(xspace.end);
    nx_.push_back(xspace.N);
    dx_.push_back(xspace.step);
    nxmin_ = std::min( nxmin_, size_t(nx_[j]) );
    nxmax_ = std::max( nxmax_, size_t(nx_[j]) );
  }
}

Structured::XSpace::Implementation::Implementation( const std::vector<Config>& config_list ) {

    reserve( config_list.size() );

    nxmin_ = std::numeric_limits<size_t>::max();
    nxmax_ = 0;

    std::string xspace_type;
    for( size_t j=0; j<ny(); ++j ) {
        config_list[j].get("type",xspace_type);
        ASSERT( xspace_type == "linear" );
        spacing::LinearSpacing::Params xspace( config_list[j] );
        xmin_.push_back(xspace.start);
        xmax_.push_back(xspace.end);
        nx_.push_back(xspace.N);
        dx_.push_back(xspace.step);
        nxmin_ = std::min( nxmin_, size_t(nx_[j]) );
        nxmax_ = std::max( nxmax_, size_t(nx_[j]) );
    }
}

void Structured::XSpace::Implementation::Implementation::reserve( long ny ) {
  ny_ = ny;
  nx_.  reserve(ny);
  xmin_.reserve(ny);
  xmax_.reserve(ny);
  dx_  .reserve(ny);
}


Structured::XSpace::Implementation::Implementation( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint ) :
  ny_(N.size()),
  nx_(N),
  xmin_(ny_,interval[0]),
  xmax_(ny_,interval[1]),
  dx_(ny_) {
  nxmin_ = std::numeric_limits<size_t>::max();
  nxmax_ = 0;
  double length = interval[1] - interval[0];
  for( size_t j=0; j<ny_; ++j ) {
    nxmin_ = std::min( nxmin_, size_t(nx_[j]) );
    nxmax_ = std::max( nxmax_, size_t(nx_[j]) );
    dx_[j] = endpoint ? length/double(nx_[j]-1) : length/double(nx_[j]);
  }
}

Structured::XSpace::Implementation::Implementation( const Spacing& spacing ) :
  ny_(1),
  nx_(ny_,spacing.size()),
  xmin_(ny_,spacing.min()),
  xmax_(ny_,spacing.max()),
  dx_(ny_) {
  const spacing::LinearSpacing& linspace = dynamic_cast<const spacing::LinearSpacing&>(*spacing.get());
  dx_[0] = linspace.step();
  nxmax_ = nx_[0];
  nxmin_ = nx_[0];
}


Grid::Spec Structured::XSpace::Implementation::spec() const {
  Grid::Spec spec;

  bool same_xmin = true;
  bool same_xmax = true;
  bool same_nx   = true;

  double xmin = xmin_[0];
  double xmax = xmax_[0];
  long   nx   = nx_  [0];
  double dx   = dx_  [0];

  ASSERT(xmin_.size() == ny_);
  ASSERT(xmax_.size() == ny_);
  ASSERT(nx_  .size() == ny_);

  for( size_t j=1; j<ny_; ++j ) {
    same_xmin = same_xmin && ( xmin_[j] == xmin );
    same_xmax = same_xmax && ( xmax_[j] == xmax );
    same_nx   = same_nx   && ( nx_  [j] == nx   );
  }

  bool endpoint = std::abs( (xmax - xmin) - (nx-1)*dx ) < 1.e-10;

  spec.set("type","linear");
  if( same_xmin ) {
    spec.set("start",xmin);
  } else {
    spec.set("start[]",xmin_);
  }
  if( same_xmax ) {
    spec.set("end",xmax);
  } else {
    spec.set("end[]",xmax_);
  }
  if( same_nx ) {
    spec.set("N",nx);
  } else {
    spec.set("N[]",nx_);
  }
  spec.set("endpoint",endpoint);

  return spec;
}

namespace {
  class Normalise {
  public:
    Normalise(const domain::RectangularDomain& domain) :
      xmin_(domain.xmin()),
      xmax_(domain.xmax()),
      degrees_(domain.units()=="degrees"),
      eps_(1e-12) {
    }

    double operator()(double x) const {
      if( degrees_ ) {
        while (eckit::types::is_strictly_greater<double>(xmin_, x, eps_)) {
          x += 360.;
        }
        while (eckit::types::is_strictly_greater<double>(x, xmax_, eps_)) {
          x -= 360.;
        }
      }
      return x;
    }

  private:
    const bool degrees_;
    const double xmin_;
    const double xmax_;
    const double eps_;
  };
}

void Structured::crop( const Domain& dom ) {

    if( dom.global() )
      return;

    ASSERT( dom.units() == projection().units() );

    auto zonal_domain = dynamic_cast<const domain::ZonalBandDomain*>(dom.get());
    auto rect_domain  = dynamic_cast<const domain::RectangularDomain*>(dom.get());

    if( zonal_domain ) {

      const double cropped_ymin = rect_domain->ymin();
      const double cropped_ymax = rect_domain->ymax();

      size_t jmin = ny();
      size_t jmax = 0;
      for( size_t j=0; j<ny(); ++j )
      {
          if( zonal_domain->contains_y(y(j)) ) {
              jmin = std::min(j, jmin);
              jmax = std::max(j, jmax);
          }
      }
      size_t cropped_ny = jmax-jmin+1;
      std::vector<double> cropped_y   ( y_   .begin()+jmin, y_   .begin()+jmin+cropped_ny );
      std::vector<double> cropped_xmin( xmin_.begin()+jmin, xmin_.begin()+jmin+cropped_ny );
      std::vector<double> cropped_xmax( xmax_.begin()+jmin, xmax_.begin()+jmin+cropped_ny );
      std::vector<double> cropped_dx  ( dx_  .begin()+jmin, dx_  .begin()+jmin+cropped_ny );
      std::vector<long>   cropped_nx  ( nx_  .begin()+jmin, nx_  .begin()+jmin+cropped_ny );
      ASSERT( cropped_nx.size() == cropped_ny );

      size_t cropped_nxmin, cropped_nxmax;
      cropped_nxmin = cropped_nxmax = static_cast<size_t>(cropped_nx.front());
      for (size_t j=1; j<cropped_ny; ++j) {
          cropped_nxmin = std::min(static_cast<size_t>(cropped_nx[j]),cropped_nxmin);
          cropped_nxmax = std::max(static_cast<size_t>(cropped_nx[j]),cropped_nxmax);
      }
      size_t cropped_npts = size_t(std::accumulate(cropped_nx.begin(), cropped_nx.end(), 0));

      Spacing cropped_yspace( new spacing::CustomSpacing(cropped_ny, cropped_y.data(), {cropped_ymin, cropped_ymax}) );

      // Modify grid
      {
        domain_ = dom;
        yspace_ = cropped_yspace;
        xmin_   = cropped_xmin;
        xmax_   = cropped_xmax;
        dx_     = cropped_dx;
        nx_     = cropped_nx;
        nxmin_  = cropped_nxmin;
        nxmax_  = cropped_nxmax;
        npts_   = cropped_npts;
        y_      = cropped_y;
      }

    } else if ( rect_domain ) {

        const double cropped_ymin = rect_domain->ymin();
        const double cropped_ymax = rect_domain->ymax();

        // Cropping in Y
        size_t jmin = ny();
        size_t jmax = 0;
        for( size_t j=0; j<ny(); ++j ) {
            if( rect_domain->contains_y(y(j)) ) {
                jmin = std::min(j, jmin);
                jmax = std::max(j, jmax);
            }
        }
        size_t cropped_ny = jmax-jmin+1;
        std::vector<double> cropped_y   ( y_ .begin()+jmin, y_ .begin()+jmin+cropped_ny );
        std::vector<double> cropped_dx  ( dx_.begin()+jmin, dx_.begin()+jmin+cropped_ny );

        std::vector<double> cropped_xmin( cropped_ny,  std::numeric_limits<double>::max() );
        std::vector<double> cropped_xmax( cropped_ny, -std::numeric_limits<double>::max() );
        std::vector<long>   cropped_nx  ( cropped_ny );

        // Cropping in X
        Normalise normalise(*rect_domain);
        for( size_t j=jmin, jcropped=0; j<=jmax; ++j, ++jcropped ) {
            size_t n=0;
            for( size_t i=0; i<nx(j); ++i ) {
                const double _x = normalise( x(i,j) );
                if( rect_domain->contains_x(_x) ) {
                    cropped_xmin[jcropped] = std::min( cropped_xmin[jcropped], _x );
                    cropped_xmax[jcropped] = std::max( cropped_xmax[jcropped], _x );
                    ++n;
               }
            }
            cropped_nx[jcropped] = n;
        }

        // Complete structures

        size_t cropped_nxmin, cropped_nxmax;
        cropped_nxmin = cropped_nxmax = static_cast<size_t>(cropped_nx.front());

        for (size_t j=1; j<cropped_ny; ++j) {
            cropped_nxmin = std::min(static_cast<size_t>(cropped_nx[j]),cropped_nxmin);
            cropped_nxmax = std::max(static_cast<size_t>(cropped_nx[j]),cropped_nxmax);
        }
        size_t cropped_npts = size_t(std::accumulate(cropped_nx.begin(), cropped_nx.end(), 0));

        Spacing cropped_yspace( new spacing::CustomSpacing(cropped_ny, cropped_y.data(), {cropped_ymin, cropped_ymax}) );

        // Modify grid
        {
          domain_ = dom;
          yspace_ = cropped_yspace;
          xmin_   = cropped_xmin;
          xmax_   = cropped_xmax;
          dx_     = cropped_dx;
          nx_     = cropped_nx;
          nxmin_  = cropped_nxmin;
          nxmax_  = cropped_nxmax;
          npts_   = cropped_npts;
          y_      = cropped_y;
        }
    }
    else {
      std::stringstream errmsg;
      errmsg << "Cannot crop the grid with domain " << dom;
      eckit::BadParameter(errmsg.str(), Here());
    }
}


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

void Structured::computeTruePeriodicity() {

  if( projection_.strictlyRegional() ) {
    periodic_x_ = false;

  } else if( domain_.global() ) {
    periodic_x_ = true;

  } else {
    // domain could be zonal band

    size_t j = ny()/2;
    if( xmin_[j] + (nx_[j]-1) * dx_[j] == xmax_[j] ) {
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
  return static_type();
}

void Structured::hash(eckit::Hash& h) const {

    h.add(y().data(), sizeof(double)*y().size());
    h.add(nx().data(), sizeof(long)*ny());

    // also add lonmin and lonmax
    h.add(xmin_.data(), sizeof(double)*xmin_.size());
    h.add(dx_.data(),   sizeof(double)*dx_.size());

    // also add projection information
    projection().hash(h);

    // also add domain information, even though already encoded in grid.
    domain().hash(h);
}

Grid::Spec Structured::spec() const {
    Grid::Spec grid_spec;

    if( name() == "structured" ) {
        grid_spec.set("type", type());
        grid_spec.set("xspace",xspace().spec());
        grid_spec.set("yspace",yspace().spec());
    }
    else {
      grid_spec.set("name",name());
    }
    grid_spec.set("domain",domain().spec());
    grid_spec.set("projection",projection().spec());
    return grid_spec;
}


// --------------------------------------------------------------------

namespace { // anonymous

static class structured : public GridBuilder {

  using Implementation = atlas::Grid::Implementation;
  using Config = Grid::Config;
  using XSpace = StructuredGrid::XSpace;

public:

    structured(): GridBuilder( "structured" ){}

    virtual void print(std::ostream& os) const {
        os << std::left << std::setw(20) << " " << "Structured grid";
    }

    virtual const Implementation* create( const std::string& name, const Config& config ) const {
        throw eckit::NotImplemented( "Cannot create structured grid from name", Here() );
    }

    virtual const Implementation* create( const Config& config ) const {

        Projection projection;
        Spacing    yspace;
        Domain     domain;

        Config config_proj;
        if( config.get("projection",config_proj) )
            projection = Projection(config_proj);

        Config config_domain;
        if( config.get("domain",config_domain) ) {
            domain = Domain(config_domain);
        }

        Config config_yspace;
        if( not config.get("yspace",config_yspace) )
            throw eckit::BadParameter("yspace missing in configuration");
        yspace = Spacing(config_yspace);

        size_t ny = yspace.size();

        XSpace xspace;

        Config config_xspace;
        std::vector<Config> config_xspace_list;

        if( config.get("xspace[]",config_xspace_list) ) {
            xspace = XSpace( config_xspace_list );
        } else if( config.get("xspace",config_xspace) ) {
            xspace = XSpace( config_xspace );
        } else {
            throw eckit::BadParameter("xspace missing in configuration");
        }

        return new StructuredGrid::grid_t(xspace, yspace, projection, domain );
    }

} structured_;

} // anonymous namespace


// --------------------------------------------------------------------


extern "C" {


    long atlas__grid__Structured__ny(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->ny();
        );
        return 0;
    }


    long atlas__grid__Structured__nx(Structured* This, long jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nx(jlat);
        );
        return 0;
    }


    void atlas__grid__Structured__nx_array(Structured* This, const long* &nx_array, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            nx_array = This->nx().data();
            size  = This->nx().size();
        );
    }


    long atlas__grid__Structured__nxmax(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nxmax();
        );
        return 0;
    }


    long atlas__grid__Structured__nxmin(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nxmin();
        );
        return 0;
    }


    long atlas__grid__Structured__size(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->size();
        );
        return 0;
    }


    double atlas__grid__Structured__y(Structured* This, long j) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->y(j);
        );
        return 0.;
    }


    double atlas__grid__Structured__x( Structured* This, long i, long j ) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->x(i,j);
        );
        return 0.;
    }


    void atlas__grid__Structured__xy( Structured* This, long i, long j, double crd[] ) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            This->xy(i,j, crd);
        );
    }

    void atlas__grid__Structured__lonlat( Structured* This, long i, long j, double crd[] ) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            This->lonlat(i,j, crd);
        );
    }


    void atlas__grid__Structured__y_array( Structured* This, const double* &y_array, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            y_array  = This->y().data();
            size     = This->y().size();
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


}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

