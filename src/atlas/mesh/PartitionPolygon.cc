/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <sstream>

#include "eckit/exception/Exceptions.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh.h"
#include "atlas/mesh/PartitionPolygon.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {

namespace {
util::Polygon::edge_set_t compute_edges(const detail::MeshImpl& mesh, idx_t halo) {
    ATLAS_TRACE("PartitionPolygon [halo=" + std::to_string(halo) + "]");
    // extract partition boundary edges by always attempting first to`
    // remove a reversed edge from a neighbouring element, if any

    util::Polygon::edge_set_t edges;
    for (idx_t t = 0; t < mesh.cells().nb_types(); ++t) {
        const Elements& elements = mesh.cells().elements(t);

        const BlockConnectivity& conn = elements.node_connectivity();
        auto field_flags              = elements.view<int, 1>(elements.flags());
        auto field_halo               = elements.view<int, 1>(elements.halo());

        constexpr bool include_patch = false;
        auto include                 = [&](idx_t e) {
            using util::Topology;
            auto topology = Topology::view(field_flags(e));
            if (field_halo(e) > halo) {
                return false;
            }
            if (not include_patch) {
                if (topology.check(Topology::PATCH)) {
                    return false;
                }
            }
            return true;
        };


        const idx_t nb_nodes = elements.nb_nodes();

        for (idx_t j = 0; j < elements.size(); ++j) {
            if (include(j)) {
                for (idx_t k = 0; k < nb_nodes; ++k) {
                    util::Polygon::edge_t edge(conn(j, k), conn(j, (k + 1) % nb_nodes));
                    if (!edges.erase(edge.reverse())) {
                        edges.insert(edge);
                    }
                }
            }
        }
    }
    return edges;
}
}  // namespace

PartitionPolygon::PartitionPolygon(const detail::MeshImpl& mesh, idx_t halo):
    util::PartitionPolygon(compute_edges(mesh, halo)), mesh_(mesh), halo_(halo) {}

size_t PartitionPolygon::footprint() const {
    size_t size = sizeof(*this);
    size += capacity() * sizeof(idx_t);
    return size;
}

void PartitionPolygon::outputPythonScript(const eckit::PathName& filepath, const eckit::Configuration& config) const {
    const auto& comm = mpi::comm(mesh_.mpi_comm());
    int mpi_rank     = int(comm.rank());
    int mpi_size     = int(comm.size());

    std::string coordinates = config.getString("coordinates", "xy");
    if (coordinates == "ij") {
        coordinates = "xy";
    }
    auto points =
        coordinates == "xy" ? this->xy() :
                            // coordinates == "ij" ? this->ij() : /* no member named `ij` in PartitionPolygon */
            this->lonlat();

    ATLAS_ASSERT(points.size() == size());
    const auto nodes_xy = array::make_view<double, 2>(mesh_.nodes().field(coordinates));
    const auto nodes_g  = array::make_view<gidx_t, 1>(mesh_.nodes().global_index());
    double xmin         = std::numeric_limits<double>::max();
    double xmax         = -std::numeric_limits<double>::max();
    for (const auto& p : points) {
        xmin = std::min(xmin, p[XX]);
        xmax = std::max(xmax, p[XX]);
    }
    comm.allReduceInPlace(xmin, eckit::mpi::min());
    comm.allReduceInPlace(xmax, eckit::mpi::max());

    idx_t count     = mesh_.nodes().size();
    idx_t count_all = count;
    comm.allReduceInPlace(count_all, eckit::mpi::sum());

    bool plot_nodes = config.getBool("nodes", false);
    for (int r = 0; r < mpi_size; ++r) {
        if (mpi_rank == r) {
            std::ofstream f(filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app);
            if (!f.is_open()) {
                throw eckit::CantOpenFile(filepath);
            }
            //clang-format off
            if (mpi_rank == 0) {
                f << "\n"
                     "# Configuration option to plot nodes"
                     "\n"
                     "plot_nodes = " +
                         std::string(plot_nodes ? "True" : "False")
                  << "\n"
                     "plot_node_numbers = False"
                     "\n"
                     "\n"
                     "import matplotlib.pyplot as plt"
                     "\n"
                     "from matplotlib.path import Path"
                     "\n"
                     "import matplotlib.patches as patches"
                     "\n"
                     "\n"
                     "from itertools import cycle"
                     "\n"
                     "import matplotlib.cm as cm"
                     "\n"
                     "import numpy as np"
                     "\n"
                     "colours = cycle([cm.Paired(i) for i in np.linspace(0,1,12,endpoint=True)])"
                     "\n"
                     "\n"
                     "fig = plt.figure()"
                     "\n"
                     "ax = fig.add_subplot(111,aspect='equal')"
                     "\n";
            }
            f << "\n"
                 "verts_"
              << r << " = [";
            for (const auto& p : points) {
                f << "\n  (" << p[XX] << ", " << p[YY] << "), ";
            }
            f << "\n]"
                 "\n"
                 "\n"
                 "codes_"
              << r
              << " = [Path.MOVETO]"
                 "\n"
                 "codes_"
              << r << ".extend([Path.LINETO] * " << (size() - 2)
              << ")"
                 "\n"
                 "codes_"
              << r
              << ".extend([Path.CLOSEPOLY])"
                 "\n"
                 "\n"
                 "count_"
              << r << " = " << count
              << "\n"
                 "count_all_"
              << r << " = " << count_all << "\n";
            if (plot_nodes) {
                f << "\n"
                     "g_"
                  << r << " = [";
                for (idx_t i = 0; i < count; ++i) {
                    f << nodes_g(i) << ", ";
                }
                f << "]"
                     "\n"
                     "x_"
                  << r << " = [";
                for (idx_t i = 0; i < count; ++i) {
                    f << nodes_xy(i, XX) << ", ";
                }
                f << "]"
                     "\n"
                     "y_"
                  << r << " = [";
                for (idx_t i = 0; i < count; ++i) {
                    f << nodes_xy(i, YY) << ", ";
                }
                f << "]";
            }
            f << "\n"
                 "\n"
                 "c = next(colours)"
                 "\n"
                 "ax.add_patch(patches.PathPatch(Path(verts_"
              << r << ", codes_" << r << "), color=c, alpha=0.3, lw=1))";
            if (plot_nodes) {
                f << "\n"
                     "if plot_nodes:"
                     "\n"
                     "    ax.scatter(x_"
                  << r << ", y_" << r
                  << ", color=c, marker='o')"
                     "\n"
                     "    if plot_node_numbers:"
                     "\n"
                     "        for (g,x,y) in zip(g_"
                  << r << ", x_" << r << ", y_" << r
                  << "):"
                     "\n"
                     "            ax.text(x,y, str(g), color=c, fontsize=12)";
            }
            f << "\n";
            // clang-format on
            if (mpi_rank == mpi_size - 1) {
                auto xticks = [&]() -> std::string {
                    std::vector<int> xticks;
                    xticks.push_back(std::floor(xmin - int(xmin) % 45));
                    while (xticks.back() < int(std::round(xmax))) {
                        xticks.push_back(xticks.back() + 45);
                    }
                    std::stringstream xticks_ss;
                    for (int i = 0; i < xticks.size(); ++i) {
                        xticks_ss << (i == 0 ? "[" : ",") << xticks[i];
                    }
                    xticks_ss << "]";
                    return xticks_ss.str();
                };
                // clang-format off
                f << "\n" "ax.set_xlim( " << xmin << "-5, " << xmax << "+5)"
                     "\n" "ax.set_ylim(-90-5, 90+5)"
                     "\n" "ax.set_xticks(" << xticks() << ")"
                     "\n" "ax.set_yticks([-90,-45,0,45,90])"
                     "\n" "plt.grid()"
                     "\n" "plt.show()";
                // clang-format on
            }
        }
        comm.barrier();
    }
}

void PartitionPolygon::allGather(util::PartitionPolygons& polygons) const {
    ATLAS_TRACE();

    const auto& comm   = mpi::comm(mesh_.mpi_comm());
    const int mpi_size = int(comm.size());

    polygons.clear();
    polygons.reserve(mpi_size);

    auto& poly = *this;

    std::vector<double> mypolygon;
    mypolygon.reserve(poly.size() * 2);

    auto points_xy = poly.xy();
    ATLAS_ASSERT(points_xy.front() == points_xy.back());
    for (auto& p : points_xy) {
        mypolygon.push_back(p[XX]);
        mypolygon.push_back(p[YY]);
    }
    ATLAS_ASSERT(mypolygon.size() >= 4);

    eckit::mpi::Buffer<double> recv_polygons(mpi_size);

    comm.allGatherv(mypolygon.begin(), mypolygon.end(), recv_polygons);

    for (idx_t p = 0; p < mpi_size; ++p) {
        PointsXY recv_points;
        recv_points.reserve(recv_polygons.counts[p]);
        for (idx_t j = 0; j < recv_polygons.counts[p] / 2; ++j) {
            PointXY pxy(*(recv_polygons.begin() + recv_polygons.displs[p] + 2 * j + XX),
                        *(recv_polygons.begin() + recv_polygons.displs[p] + 2 * j + YY));
            recv_points.push_back(pxy);
        }
        polygons.emplace_back(new util::ExplicitPartitionPolygon(std::move(recv_points)));
    }
}

void PartitionPolygon::print(std::ostream& out) const {
    out << "polygon:{"
        << "halo:" << halo_ << ",size:" << size() << ",nodes:" << static_cast<const util::Polygon&>(*this) << "}";
}

namespace {
PartitionPolygon::PointsXY points(const Mesh::Implementation& mesh_, const Field& coords,
                                  const std::vector<idx_t>& polygon) {
    PartitionPolygon::PointsXY points_xy;
    points_xy.reserve(polygon.size());

    auto xy_view = array::make_view<double, 2>(coords);
    auto flags   = array::make_view<int, 1>(mesh_.nodes().flags());

    bool domain_includes_north_pole = false;
    bool domain_includes_south_pole = false;
    if (mesh_.grid()) {
        if (mesh_.grid().domain()) {
            domain_includes_north_pole = mesh_.grid().domain().containsNorthPole();
            domain_includes_south_pole = mesh_.grid().domain().containsSouthPole();
        }
    }
    auto bc_north = [&flags, &domain_includes_north_pole](idx_t n) {
        if (domain_includes_north_pole) {
            using Topology = atlas::mesh::Nodes::Topology;
            return Topology::check(flags(n), Topology::BC | Topology::NORTH);
        }
        return false;
    };
    auto bc_south = [&flags, &domain_includes_south_pole](idx_t n) {
        if (domain_includes_south_pole) {
            using Topology = atlas::mesh::Nodes::Topology;
            return Topology::check(flags(n), Topology::BC | Topology::SOUTH);
        }
        return false;
    };

    if (coords.name() == "xy") {
        for (idx_t i : polygon) {
            double y = bc_north(i) ? 90. : bc_south(i) ? -90. : xy_view(i, idx_t(YY));
            points_xy.emplace_back(xy_view(i, idx_t(XX)), y);
        }
    }
    else {
        for (idx_t i : polygon) {
            points_xy.emplace_back(xy_view(i, idx_t(XX)), xy_view(i, idx_t(YY)));
        }
    }
    return points_xy;
}
}  // namespace

PartitionPolygon::PointsXY PartitionPolygon::xy() const {
    return points(mesh_, mesh_.nodes().xy(), *this);
}

PartitionPolygon::PointsLonLat PartitionPolygon::lonlat() const {
    return points(mesh_, mesh_.nodes().lonlat(), *this);
}

}  // namespace mesh
}  // namespace atlas
