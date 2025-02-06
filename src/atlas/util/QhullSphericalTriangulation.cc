#include "QhullSphericalTriangulation.h"

#include "atlas/library/defines.h"

#include <cmath>
#include <random>
#include <algorithm>
#include <utility>

#if ATLAS_HAVE_QHULL

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wtemplate-id-cdtor"
#endif

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/QhullPoints.h>

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif

#include "atlas/runtime/Exception.h"

namespace atlas{
namespace util{

#if ATLAS_HAVE_QHULL
struct QhullSphericalTriangulation::Qhull : public orgQhull::Qhull {
    using orgQhull::Qhull::Qhull;
};
#else
struct QhullSphericalTriangulation::Qhull {
    template<typename... Args>
    Qhull(Args...) {
        throw_Exception("Atlas has not been compiled with Qhull",Here());
    }
};
#endif

QhullSphericalTriangulation::~QhullSphericalTriangulation() = default;

QhullSphericalTriangulation::QhullSphericalTriangulation(size_t N, const double lonlat[]) :
  QhullSphericalTriangulation(N, lonlat, lonlat+1, 2, 2) {
}

QhullSphericalTriangulation::QhullSphericalTriangulation(size_t N, const double lon[], const double lat[]) :
  QhullSphericalTriangulation(N, lon, lat, 1, 1) {
}

QhullSphericalTriangulation::QhullSphericalTriangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride) {
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

    constexpr int dim = 3;
    constexpr const char* command = "Qt";
    constexpr const char* comment = "";
    const double* coordinates = reinterpret_cast<const double*>(points_xyz_.data());
    qhull_ = std::make_unique<Qhull>(comment, dim, N, coordinates, command);
}

size_t QhullSphericalTriangulation::size() const {
#if ATLAS_HAVE_QHULL
    return qhull_->facetList().size();
#else
    return 0;
#endif
}


template <typename Qhull, typename Points, typename Value>
static inline void qhull_get_triangles(const Qhull& qhull, const Points& points_xyz_, std::array<Value,3> triangles[]) {
#if ATLAS_HAVE_QHULL
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
    for (const auto& facet : qhull.facetList()){
        auto& tri = triangles[jtri++];
        size_t jvrt{0};
        for( const auto& vertex: facet.vertices()){
            tri[jvrt++] = vertex.point().id();
        }

        ensure_outward_normal(tri);
    }
#endif
}

void QhullSphericalTriangulation::triangles(std::array<int,3> triangles[]) const {
    qhull_get_triangles(*qhull_,points_xyz_,triangles);
}

void QhullSphericalTriangulation::triangles(std::array<long,3> triangles[]) const {
    qhull_get_triangles(*qhull_,points_xyz_,triangles);
}

void QhullSphericalTriangulation::triangles(std::array<long long,3> triangles[]) const {
    qhull_get_triangles(*qhull_,points_xyz_,triangles);
}

void QhullSphericalTriangulation::triangles(std::array<size_t,3> triangles[]) const {
    qhull_get_triangles(*qhull_,points_xyz_,triangles);
}

std::vector<std::array<idx_t,3>> QhullSphericalTriangulation::triangles() const {
    return triangles<idx_t>();
}

} // namespace util
} // namespace atlas
