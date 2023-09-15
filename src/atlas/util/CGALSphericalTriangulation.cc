#include "CGALSphericalTriangulation.h"

#include "atlas/library/defines.h"

#include <cmath>
#include <random>
#include <algorithm>
#include <utility>

#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_CGAL

// CGAL needs -DCGAL_NDEBUG to reach peak performance ...
#define CGAL_NDEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>


using K            = ::CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron_3 = ::CGAL::Polyhedron_3<K>;
using Point_3      = Polyhedron_3::Point_3;

#endif


namespace atlas{
namespace util{

#if ATLAS_HAVE_CGAL
struct CGALSphericalTriangulation::CGAL {
    CGAL(const std::vector<std::array<double,3>>& pts) {
        std::vector<Point_3> vertices(pts.size());
        for (size_t i = 0, size = vertices.size(); i < size; ++i) {
            vertices[i] = Point_3(pts[i][0], pts[i][1], pts[i][2]);
            point_map_[vertices[i]] = i;
        }

        // compute convex hull of non-collinear points
        ::CGAL::convex_hull_3(vertices.begin(), vertices.end(), poly_);
    }
    std::unordered_map<Point_3,size_t> point_map_;
    Polyhedron_3 poly_;
};
#else
struct CGALSphericalTriangulation::CGAL {
    template<typename... Args>
    CGAL(Args...) {
        throw_Exception("Atlas has not been compiled with CGAL",Here());
    }
};
#endif

CGALSphericalTriangulation::~CGALSphericalTriangulation() = default;

CGALSphericalTriangulation::CGALSphericalTriangulation(size_t N, const double lonlat[]) :
  CGALSphericalTriangulation(N, lonlat, lonlat+1, 2, 2) {
}

CGALSphericalTriangulation::CGALSphericalTriangulation(size_t N, const double lon[], const double lat[]) :
  CGALSphericalTriangulation(N, lon, lat, 1, 1) {
}

CGALSphericalTriangulation::CGALSphericalTriangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride) {
    auto lonlat2xyz = [](double lon, double lat, auto& xyz) {
        constexpr double deg2rad = M_PI / 180.;
        const double lambda     = deg2rad * lon;
        const double phi        = deg2rad * lat;

        const double sin_phi    = std::sin(phi);
        const double cos_phi    = std::cos(phi);
        const double sin_lambda = std::sin(lambda);
        const double cos_lambda = std::cos(lambda);

        xyz[0] = cos_phi * cos_lambda;
        xyz[1] = cos_phi * sin_lambda;
        xyz[2] = sin_phi;
    };

    points_xyz_.resize(N);

    for (size_t i = 0; i < N; ++i) {
        lonlat2xyz(lon[i*lon_stride], lat[i*lat_stride], points_xyz_[i]);
    }

    cgal_ = std::make_unique<CGAL>(points_xyz_);
}

size_t CGALSphericalTriangulation::size() const {
#if ATLAS_HAVE_CGAL
    return cgal_->poly_.size_of_facets();
#else
    return 0;
#endif
}


template <typename CGAL, typename Points, typename Value>
static inline void CGAL_get_triangles(const CGAL& cgal, const Points& points_xyz_, std::array<Value,3> triangles[]) {
#if ATLAS_HAVE_CGAL
    auto ensure_outward_normal = [&](auto& tri) {

        auto dot = [](const auto& p1, const auto& p2) {
            return p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
        };
        auto cross = [](const auto& p1, const auto& p2) {
            return std::array<double,3>{
                p1[1] * p2[2] - p1[2] * p2[1], p1[2] * p2[0] - p1[0] * p2[2],
                p1[0] * p2[1] - p1[1] * p2[0]
            };
        };

        const std::array<double,3>& a = points_xyz_[tri[0]];
        const std::array<double,3>& b = points_xyz_[tri[1]];
        const std::array<double,3>& c = points_xyz_[tri[2]];
        std::array<double,3> ba {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
        std::array<double,3> bc {c[0]-b[0], c[1]-b[1], c[2]-b[2]};

        bool outward = dot(b, cross(bc,ba)) > 0;

        if (not outward) {
            std::swap(tri[1], tri[2]);
        }
    };

    size_t jtri{0};

    const auto& poly = cgal.poly_;
    const auto& point_map = cgal.point_map_;
    const idx_t nb_points = idx_t(points_xyz_.size());

    /* triangles */

    const idx_t nb_triags = poly.size_of_facets();

    idx_t idx[3];
    Polyhedron_3::Vertex_const_handle vts[3];

    idx_t tidx = 0;
    for (Polyhedron_3::Facet_const_iterator f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        // loop  over half-edges and take each vertex()

        auto& tri = triangles[jtri++];

        size_t jvrt{0};

        idx_t iedge                                               = 0;
        Polyhedron_3::Halfedge_around_facet_const_circulator edge = f->facet_begin();
        do {
            Polyhedron_3::Vertex_const_handle vh = edge->vertex();
            const Polyhedron_3::Point_3& p       = vh->point();
                        tri[jvrt++] = point_map.at(vh->point());

            ++iedge;
            ++edge;
        } while (edge != f->facet_begin() && iedge < 3);

        ATLAS_ASSERT(iedge == 3);

        ensure_outward_normal(tri);
    }

    ATLAS_ASSERT(jtri == nb_triags);
#endif
}

void CGALSphericalTriangulation::triangles(std::array<int,3> triangles[]) const {
    CGAL_get_triangles(*cgal_,points_xyz_,triangles);
}

void CGALSphericalTriangulation::triangles(std::array<long,3> triangles[]) const {
    CGAL_get_triangles(*cgal_,points_xyz_,triangles);
}

void CGALSphericalTriangulation::triangles(std::array<long long,3> triangles[]) const {
    CGAL_get_triangles(*cgal_,points_xyz_,triangles);
}

void CGALSphericalTriangulation::triangles(std::array<size_t,3> triangles[]) const {
    CGAL_get_triangles(*cgal_,points_xyz_,triangles);
}

std::vector<std::array<idx_t,3>> CGALSphericalTriangulation::triangles() const {
    return triangles<idx_t>();
}

} // namespace util
} // namespace atlas
