/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/FiniteElement.h"

#include "eckit/log/Plural.h"
#include "eckit/log/Seconds.h"
#include "eckit/log/Timer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/geometry/Point3.h"

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/LibAtlas.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {


MethodBuilder<FiniteElement> __builder("finite-element");


// epsilon used to scale edge tolerance when projecting ray to intesect element
static const double parametricEpsilon = 1e-16;


}  // (anonymous namespace)


void FiniteElement::setup(mesh::Mesh& meshSource, mesh::Mesh& meshTarget) {
    using namespace atlas;
    eckit::TraceTimer<LibAtlas> tim("atlas::interpolation::method::FiniteElement::setup()");


    // generate 3D point coordinates
    mesh::actions::BuildXYZField("xyz")(meshSource);
    mesh::actions::BuildXYZField("xyz")(meshTarget);


    // generate barycenters of each triangle & insert them on a kd-tree
    mesh::actions::BuildCellCentres()(meshSource);

    eckit::ScopedPtr<ElemIndex3> eTree(create_element_centre_index(meshSource));


    const mesh::Nodes  &i_nodes  = meshSource.nodes();
    icoords_ = array::ArrayView<double,2>(i_nodes.field( "xyz" ));
    ilonlat_ = array::ArrayView<double,2>(i_nodes.lonlat());

    ocoords_ = array::ArrayView<double,2>(meshTarget.nodes().field("xyz"));
    olonlat_ = array::ArrayView<double,2>(meshTarget.nodes().lonlat());

    connectivity_ = &meshSource.cells().node_connectivity();

    size_t inp_npts  = i_nodes.size();
    size_t out_npts  = meshTarget.nodes().size();

    array::ArrayView<int, 1> out_ghosts(meshTarget.nodes().ghost());

    size_t Nelements = meshSource.cells().size();
    const double maxFractionElemsToTry = 0.2;


    // weights -- one per vertex of element, triangles (3) or quads (4)

    std::vector< eckit::linalg::Triplet > weights_triplets;  // structure to fill-in sparse matrix
    weights_triplets.reserve( out_npts * 4 );                // preallocate space as if all elements where quads

    // search nearest k cell centres

    const size_t maxNbElemsToTry = 20;//std::max<size_t>(64, size_t(Nelements * maxFractionElemsToTry));
    size_t max_neighbours = 0;

    std::vector<size_t> failures;

    {
        eckit::TraceTimer<LibAtlas> timerProj("Projecting");

        for ( size_t ip = 0; ip < out_npts; ++ip ) {
            if (out_ghosts(ip)) {
                continue;
            }

            if (ip && (ip % 1000 == 0)) {
                double rate = ip / timerProj.elapsed();
                Log::debug<LibAtlas>() << eckit::BigNum(ip) << " (at " << rate << " points/s)..." << std::endl;
            }

            Point p ( ocoords_[ip].data() ); // lookup point

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
                Log::debug<LibAtlas>() << "---------------------------------------------------------------------------\n";
                Log::debug<LibAtlas>() << "Failed to project point (lon,lat)=("<<olonlat_(ip,LON)<<","<<olonlat_(ip,LAT) << ")" << '\n';
                Log::debug<LibAtlas>() << failures_log.str();
            }
        }
    }

    Log::debug<LibAtlas>() << "Projected " << eckit::Plural(out_npts, "point") << std::endl;
    Log::debug<LibAtlas>() << "Maximum neighbours searched was " << eckit::Plural(max_neighbours, "element") << std::endl;


    eckit::mpi::comm().barrier();
    if (failures.size()) {

        // If this fails, consider lowering atlas::grid::parametricEpsilon
        std::ostringstream msg;
        msg << "Rank " << eckit::mpi::comm().rank() << " failed to project points:\n";
        for (std::vector<size_t>::const_iterator i = failures.begin(); i != failures.end(); ++i) {
            msg << "\t(lon,lat) = (" << olonlat_[*i][LON] << "," << olonlat_[*i][LAT] << ")\n";
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

    const size_t inp_points = icoords_.shape(0);
    size_t idx[4];
    double w[4];

    Triplets triplets;
    Ray ray( ocoords_[ip].data() );

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
                    icoords_[idx[0]].data(),
                    icoords_[idx[1]].data(),
                    icoords_[idx[2]].data());

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
            else { // fallback to 2d lonlat elements

                element::Triag2D triag2d(ilonlat_[idx[0]].data(),
                                         ilonlat_[idx[1]].data(),
                                         ilonlat_[idx[2]].data());

                is = triag2d.intersects(Vector2D::Map(olonlat_[ip].data()),edgeEpsilon);

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

                if(Log::debug<LibAtlas>()) {
                    failures_log << "Failed to project point " << Point(ocoords_[ip].data()) << " to " << triag2d << " with " << is << '\n';
                }
            }


        } else {

            /* quadrilateral */
            element::Quad3D quad(
                    icoords_[idx[0]].data(),
                    icoords_[idx[1]].data(),
                    icoords_[idx[2]].data(),
                    icoords_[idx[3]].data() );

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
            else {

                if( Log::debug<LibAtlas>() ) {
                    failures_log << "Failed to project point " << Point(ocoords_[ip].data()) << " to " << quad << " with t=" << is.t << '\n';
                }
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

