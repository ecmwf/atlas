/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction. and Interpolation
 */

#include <cmath>

#include "atlas/interpolation/method/FiniteElement.h"

#include "eckit/log/Plural.h"
#include "eckit/log/Seconds.h"
#include "eckit/log/Timer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/geometry/Point3.h"

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Point.h"
#include "atlas/util/CoordinateEnums.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {


MethodBuilder<FiniteElement> __builder("finite-element");


// epsilon used to scale edge tolerance when projecting ray to intesect element
static const double parametricEpsilon = 1e-16;

inline PointLonLat to_lonlat(const PointXYZ& p) {
  static double rad2deg = 180.*M_1_PI;
  static double eps = 1.e-12;
  double z = p.z();
  if      (z >  1.0) z =  1.0;
  else if (z < -1.0) z = -1.0;
  double y = p.y();
  if( std::abs(y) < eps ) y *= 0.;
  return PointLonLat (  std::atan2(y, p.x()) *rad2deg ,
                        std::asin (z)        *rad2deg );
}

}  // (anonymous namespace)

void FiniteElement::setup(const FunctionSpace& source, const FunctionSpace& target) {
    eckit::TraceTimer<Atlas> tim("atlas::interpolation::method::FiniteElement::setup()");

    if( functionspace::NodeColumns tgt = target ) {

        Mesh meshTarget = tgt.mesh();

        // generate 3D point coordinates
        target_xyz_ = mesh::actions::BuildXYZField("xyz")(meshTarget);
        target_ghost_ = meshTarget.nodes().ghost();

    }
    else if ( functionspace::PointCloud tgt = target ) {

        const size_t N = tgt.size();
        target_xyz_   = Field("xyz",   array::make_datatype<double>(), array::make_shape(N,3) );
        target_ghost_ = tgt.ghost();
        array::ArrayView<double,2> lonlat = array::make_view<double,2>( tgt.lonlat() );
        array::ArrayView<double,2> xyz    = array::make_view<double,2>( target_xyz_ );
        PointXYZ p2;
        for( size_t n=0; n<N; ++n ) {
            const PointLonLat p1(lonlat(n, 0), lonlat(n, 1));
            util::Earth::convertGeodeticToGeocentric(p1, p2);
            xyz(n, 0) = p2.x();
            xyz(n, 1) = p2.y();
            xyz(n, 2) = p2.z();
        }

    }
    else {
        NOTIMP;
    }

    setup(source);
}
    
void FiniteElement::setup(const FunctionSpace& source) {

    const functionspace::NodeColumns src = source;
    ASSERT(src);

    Mesh meshSource = src.mesh();

    // generate 3D point coordinates
    Field source_xyz = mesh::actions::BuildXYZField("xyz")(meshSource);

    // generate barycenters of each triangle & insert them on a kd-tree
    Field cell_centres = mesh::actions::BuildCellCentres("centre")(meshSource);

    eckit::ScopedPtr<ElemIndex3> eTree( create_element_kdtree( cell_centres ) );

    const mesh::Nodes  &i_nodes  = meshSource.nodes();

    icoords_.reset( new array::ArrayView<double,2>( array::make_view<double,2>(source_xyz)) );
    ocoords_.reset( new array::ArrayView<double,2>( array::make_view<double,2>(target_xyz_)) );

    connectivity_ = &meshSource.cells().node_connectivity();

    size_t inp_npts  = i_nodes.size();
    size_t out_npts  = ocoords_->shape(0);

    array::ArrayView<int, 1> out_ghosts = array::make_view<int,1>(target_ghost_);

    size_t Nelements = meshSource.cells().size();
    const double maxFractionElemsToTry = 0.2;


    // weights -- one per vertex of element, triangles (3) or quads (4)

    std::vector< eckit::linalg::Triplet > weights_triplets;  // structure to fill-in sparse matrix
    weights_triplets.reserve( out_npts * 4 );                // preallocate space as if all elements where quads

    // search nearest k cell centres

    const size_t maxNbElemsToTry = std::max<size_t>(64, size_t(Nelements * maxFractionElemsToTry));
    size_t max_neighbours = 0;

    std::vector<size_t> failures;

    {
        std::stringstream msg; msg << "Computing interpolation weights for " << eckit::BigNum(out_npts) << " points...";
        Log::debug<Atlas>() << msg.str() << std::endl;
        eckit::TraceTimer<Atlas> timerProj(msg.str()+"done");

        for ( size_t ip = 0; ip < out_npts; ++ip ) {
            if (out_ghosts(ip)) {
                continue;
            }

            if (ip && (ip % 1000 == 0)) {
                double rate = ip / timerProj.elapsed();
                Log::debug<Atlas>() << "    " << eckit::BigNum(ip)<<" points completed (at " << rate << " points/s)..." << std::endl;
            }

            PointXYZ p ( (*ocoords_)[ip].data() ); // lookup point

            size_t kpts = 1;
            bool success = false;
            std::ostringstream failures_log;

            while (!success && kpts <= maxNbElemsToTry) {

                max_neighbours = std::max(kpts, max_neighbours);

                ElemIndex3::NodeList cs = eTree->kNearestNeighbours(p, kpts);
                Triplets triplets = projectPointToElements(ip,cs,failures_log);

                if (triplets.size()) {
                    std::copy(triplets.begin(), triplets.end(), std::back_inserter(weights_triplets));
                    success = true;
                }
                kpts *= 2;
            }

            if (!success) {
                failures.push_back(ip);
                Log::debug<Atlas>() << "---------------------------------------------------------------------------\n";
                PointLonLat pll = to_lonlat( p );
                Log::debug<Atlas>() << "Failed to project point (lon,lat)="<< pll << '\n';
                Log::debug<Atlas>() << failures_log.str();
            }
        }
    }
    Log::debug<Atlas>() << "Maximum neighbours searched was " << eckit::Plural(max_neighbours, "element") << std::endl;


    eckit::mpi::comm().barrier();
    if (failures.size()) {

        // If this fails, consider lowering atlas::grid::parametricEpsilon
        std::ostringstream msg;
        msg << "Rank " << eckit::mpi::comm().rank() << " failed to project points:\n";
        for (std::vector<size_t>::const_iterator i = failures.begin(); i != failures.end(); ++i) {
            PointXYZ p ( (*ocoords_)[*i].data() ); // lookup point
            PointLonLat pll = to_lonlat( p );
            msg << "\t(lon,lat) = " << pll << "\n";
        }


        Log::error() << msg.str() << std::endl;
        throw eckit::SeriousBug(msg.str());
    }

    // fill sparse matrix and return
    Matrix A(out_npts, inp_npts, weights_triplets);
    matrix_.swap(A);
}


Method::Triplets FiniteElement::projectPointToElements(
        size_t ip,
        const ElemIndex3::NodeList& elems,
        std::ostream& failures_log ) const {

    ASSERT(elems.begin() != elems.end());

    const size_t inp_points = icoords_->shape(0);
    size_t idx[4];
    double w[4];

    Triplets triplets;
    Ray ray( (*ocoords_)[ip].data() );

    for (ElemIndex3::NodeList::const_iterator itc = elems.begin(); itc != elems.end(); ++itc) {

        const size_t elem_id = (*itc).value().payload();
        ASSERT(elem_id < connectivity_->rows());

        const size_t nb_cols = connectivity_->cols(elem_id);
        ASSERT(nb_cols == 3 || nb_cols == 4);

        for (size_t i = 0; i < nb_cols; ++i) {
            idx[i] = size_t((*connectivity_)(elem_id, i));
            ASSERT(idx[i] < inp_points);
        }

        if (nb_cols == 3) {

            /* triangle */
            element::Triag3D triag(
                    (*icoords_)[idx[0]].data(),
                    (*icoords_)[idx[1]].data(),
                    (*icoords_)[idx[2]].data());

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(triag.area());
            ASSERT(edgeEpsilon >= 0);

            Intersect is = triag.intersects(ray, edgeEpsilon);

            if (is) {

                // weights are the linear Lagrange function evaluated at u,v (aka barycentric coordinates)
                w[0] = 1. - is.u - is.v;
                w[1] = is.u;
                w[2] = is.v;

                for (size_t i = 0; i < 3; ++i) {
                    triplets.push_back( Triplet( ip, idx[i], w[i] ) );
                }

                break; // stop looking for elements
            }

        } else {

            /* quadrilateral */
            element::Quad3D quad(
                    (*icoords_)[idx[0]].data(),
                    (*icoords_)[idx[1]].data(),
                    (*icoords_)[idx[2]].data(),
                    (*icoords_)[idx[3]].data() );

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(quad.area());
            ASSERT(edgeEpsilon >= 0);

            Intersect is = quad.intersects(ray, edgeEpsilon);

            if (is) {

                // weights are the bilinear Lagrange function evaluated at u,v
                w[0] = (1. - is.u) * (1. - is.v);
                w[1] =       is.u  * (1. - is.v);
                w[2] =       is.u  *       is.v ;
                w[3] = (1. - is.u) *       is.v ;

                for (size_t i = 0; i < 4; ++i) {
                    triplets.push_back( Triplet( ip, idx[i], w[i] ) );
                }
                break; // stop looking for elements
            }

        }

    } // loop over nearest elements

    if (!triplets.empty()) {
        normalise(triplets);
    }
    return triplets;
}


}  // method
}  // interpolation
}  // atlas

