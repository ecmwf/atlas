/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// file deepcode ignore MissingOpenCheckOnFile: False positive

#include "atlas/grid/StructuredPartitionPolygon.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace grid {

void compute(const functionspace::FunctionSpaceImpl& _fs, idx_t _halo, std::vector<Point2>& points,
             std::vector<Point2>& bb) {
    if (not dynamic_cast<const functionspace::detail::StructuredColumns*>(&_fs)) {
        throw_Exception("Could not cast functionspace to StructuredColumns", Here());
    }
    const auto& fs  = dynamic_cast<const functionspace::detail::StructuredColumns&>(_fs);
    const auto grid = fs.grid();
    const auto dom  = RectangularDomain(grid.domain());

    if (_halo > 0 && _halo != fs.halo()) {
        throw_Exception("halo must zero or matching the one of the StructuredColumns", Here());
    }

    const int y_numbering = (grid.y().front() < grid.y().back()) ? +1 : -1;
    bool jfirst_at_pole = y_numbering * 90. + grid.y(0) == 0.;
    bool jlast_at_pole  = y_numbering * 90. - grid.y(grid.ny() - 1) == 0;

    auto compute_j = [&](const idx_t j) {
        if (j < 0) {
            return -j - 1 + jfirst_at_pole;
        }
        else if (j >= grid.ny()) {
            return 2 * grid.ny() - j - 1 - jlast_at_pole;
        }
        return j;
    };

    auto compute_y = [&](const idx_t j) {
        return fs.compute_xy(0,j).y();
    };

    auto compute_x = [&](const idx_t i, const idx_t j) {
        return fs.compute_xy(i,j).x();
    };

    auto equal = [](const double& a, const double& b) { return std::abs(a - b) < 1.e-12; };

    auto last_edge_horizontal = [&]() {
        if (points.size() < 2) {
            return false;
        }
        size_t size = points.size();
        return equal(points.at(size - 1)[YY], points.at(size - 2)[YY]);
    };
    auto last_edge_vertical = [&]() {
        if (points.size() < 2) {
            return false;
        }
        size_t size = points.size();
        return equal(points.at(size - 1)[XX], points.at(size - 2)[XX]);
    };

    PointXY p;
    idx_t c{0};

    auto add_point = [&](const Point2& point) {
        points.emplace_back(point);
        // Log::info() << "add point (" << points.size() - 1 << ")  " << point << std::endl;
    };

    auto add_vertical_edge = [&](const Point2& point) {
        if (last_edge_vertical()) {
            points.back()[YY] = point[YY];
            //            Log::info() << "mod point (" << points.size() - 1 << ")  " << point << std::endl;
        }
        else {
            add_point(point);
            c++;
        }
    };
    auto add_horizontal_edge = [&](const Point2& point) {
        if (last_edge_horizontal()) {
            points.back()[XX] = point[XX];
            //            Log::info() << "mod point (" << points.size() - 1 << ")  " << point << std::endl;
        }
        else {
            add_point(point);
            c++;
        }
    };

    double ymax, ymin, xmax, xmin;
    if (fs.j_begin() >= fs.j_end()) {
        return;
    }

    // Top
    // Top left point
    {
        idx_t j = _halo ? fs.j_begin_halo() : fs.j_begin();
        idx_t i = _halo ? fs.i_begin_halo(j) : fs.i_begin(j);
        if (j == 0 && _halo == 0) {
            p[YY] = y_numbering < 0 ? dom.ymax() : dom.ymin();
        }
        else {
            p[YY] = 0.5 * (compute_y(j - 1) + compute_y(j));
        }
        if (i == 0 && _halo == 0) {
            p[XX] = dom.xmin();
        }
        else {
            p[XX] = 0.5 * (compute_x(i - 1, j) + compute_x(i, j));
        }
        add_point(p);
    }

    // Top right point
    {
        idx_t j = _halo ? fs.j_begin_halo() :  fs.j_begin();
        idx_t i = _halo ? fs.i_end_halo(j) - 1 : fs.i_end(j) - 1;
        if (j == 0 && _halo == 0) {
            p[YY] = y_numbering < 0 ? dom.ymax() : dom.ymin();
        }
        else {
            p[YY] = 0.5 * (compute_y(j - 1) + compute_y(j));
        }
        if (i == grid.nx(compute_j(j)) - 1 && _halo == 0) {
            p[XX] = dom.xmax();
        }
        else {
            p[XX] = 0.5 * (compute_x(i, j) + compute_x(i + 1, j));
        }
        add_horizontal_edge(p);

        ymax = p[YY];
        xmax = p[XX];
    }
    // Right side
    {
        idx_t j_begin = _halo ? fs.j_begin_halo() : fs.j_begin();
        idx_t j_end = _halo ? fs.j_end_halo() : fs.j_end();
        for (idx_t j = j_begin; j < j_end - 1 ; ++j) {
            p[YY] = compute_y(j);
            idx_t i = _halo ? fs.i_end_halo(j) - 1 : fs.i_end(j) - 1;
            if (i == grid.nx(compute_j(j)) - 1 && _halo == 0) {
                p[XX] = dom.xmax();
            }
            else {
                p[XX] = 0.5 * (compute_x(i, j) + compute_x(i + 1, j));
            }
            if (p == points.back()) {
                continue;
            }
            if (not equal(p[XX], points.back()[XX])) {
                // Make a corner plus horizontal edge

                PointXY ptmp = p;

                // vertical edge
                p[XX] = points.back()[XX];
                p[YY] = 0.5 * (points.back()[YY] + p[YY]);
                add_vertical_edge(p);
                p = ptmp;

                // horizontal edge
                ptmp  = p;
                p[YY] = points.back()[YY];
                add_horizontal_edge(p);
                p = ptmp;
            }
            add_vertical_edge(p);

            xmax = std::min(xmax, p[XX]);
        }
    }
    // Bottom
    // Bottom right point(s)
    {
        idx_t j = _halo ? fs.j_end_halo() - 1 : fs.j_end() - 1;
        idx_t i = _halo ? fs.i_end_halo(j) - 1 : fs.i_end(j) - 1;

        if (j == grid.ny() - 1 && _halo == 0) {
            p[YY] = y_numbering < 0 ? dom.ymin() : dom.ymax();
        }
        else {
            p[YY] = 0.5 * (compute_y(j) + compute_y(j + 1));
        }
        if (i == grid.nx(compute_j(j)) - 1 && _halo == 0) {
            p[XX] = dom.xmax();
        }
        else {
            p[XX] = 0.5 * (compute_x(i, j) + compute_x(i + 1, j));
        }

        ymin = p[YY];

        PointXY pmin = p;

        if (not equal(p[XX], points.back()[XX])) {
            PointXY ptmp;
            ptmp  = p;
            p[XX] = points.back()[XX];
            p[YY] = 0.5 * (points.back()[YY] + compute_y(j));
            add_vertical_edge(p);

            pmin = p;
            xmax = std::min(xmax, p[XX]);

            p = ptmp;


            ptmp  = p;
            p[YY] = points.back()[YY];

            add_horizontal_edge(p);
            p = ptmp;
        }
        if (xmax - grid.xspace().dx()[compute_j(j)] < compute_x(i, j)) {
            xmax = std::min(xmax, 0.5 * (compute_x(i + 1, j) + compute_x(i, j)));
        }
        else {
            ymin = pmin[YY];
        }
        add_vertical_edge(p);
    }

    // Bottom left point
    {
        idx_t j = _halo ? fs.j_end_halo() - 1 : fs.j_end() - 1;
        idx_t i = _halo ? fs.i_begin_halo(j) : fs.i_begin(j);
        if (j == grid.ny() - 1 && _halo == 0) {
            p[YY] = y_numbering < 0 ? dom.ymin() : dom.ymax();
        }
        else {
            p[YY] = 0.5 * (compute_y(j) + compute_y(j + 1));
        }
        if (i == 0 && _halo == 0) {
            p[XX] = dom.xmin();
        }
        else {
            p[XX] = 0.5 * (compute_x(i - 1, j) + compute_x(i, j));
        }
        add_horizontal_edge(p);
        xmin = p[XX];
    }
    // Left side
    {
        idx_t j_end = _halo ? fs.j_end_halo() : fs.j_end();
        idx_t j_begin = _halo ? fs.j_begin_halo() : fs.j_begin();
        for (idx_t j = j_end - 1; j >= j_begin; --j) {
            p[YY] = compute_y(j);
            idx_t i = _halo ? fs.i_begin_halo(j) : fs.i_begin(j);
            if (i == 0 && _halo == 0) {
                p[XX] = dom.xmin();
            }
            else {
                p[XX] = 0.5 * (compute_x(i - 1, j) + compute_x(i, j));
            }

            if (j > j_begin) {
                xmin = std::max(xmin, p[XX]);
            }
            if (j == j_begin && xmin < points[0][XX]) {
                idx_t jtop = j_begin;
                idx_t itop = _halo ? fs.i_begin_halo(jtop) : fs.i_begin(jtop);
                if (xmin + grid.xspace().dx()[compute_j(jtop)] > compute_x(itop, jtop)) {
                    xmin = std::max(xmin, 0.5 * (compute_x(itop - 1, jtop) + compute_x(itop, jtop)));
                }
                else {
                    ymax = 0.5 * (p[YY] + compute_y(jtop+1));
                }
            }

            if (p == points.back()) {
                continue;
            }

            if (not equal(p[XX], points.back()[XX])) {
                PointXY ptmp;
                ptmp = p;

                // vertical edge
                p[XX] = points.back()[XX];
                p[YY] = 0.5 * (points.back()[YY] + compute_y(j));
                add_vertical_edge(p);
                p = ptmp;

                // horizontal edge
                ptmp  = p;
                p[YY] = points.back()[YY];
                add_horizontal_edge(p);
                p = ptmp;
            }

            add_vertical_edge(p);
        }
    }
    // Connect to top
    if (last_edge_vertical()) {
        points.pop_back();
        c--;
    }
    add_point(points.front());

    bb = std::vector<Point2>{{xmin, ymax}, {xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
}

StructuredPartitionPolygon::StructuredPartitionPolygon(const functionspace::FunctionSpaceImpl& fs, idx_t halo):
    fs_(fs), halo_(halo) {
    ATLAS_TRACE("StructuredPartitionPolygon");
    compute(fs, halo, points_, inner_bounding_box_);
    auto min = Point2(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    auto max = Point2(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    for (int i = 0; i < static_cast<int>(inner_bounding_box_.size()) - 1; ++i) {
        min = Point2::componentsMin(min, inner_bounding_box_[i]);
        max = Point2::componentsMax(max, inner_bounding_box_[i]);
    }
    inscribed_domain_ = {{min[XX], max[XX]}, {min[YY], max[YY]}};
}

size_t StructuredPartitionPolygon::footprint() const {
    size_t size = sizeof(*this);
    size += capacity() * sizeof(idx_t);
    return size;
}

void StructuredPartitionPolygon::outputPythonScript(const eckit::PathName& filepath,
                                                    const eckit::Configuration& config) const {
    ATLAS_TRACE("Output PartitionPolygon");
    const auto& comm = atlas::mpi::comm(fs_.mpi_comm());
    int mpi_rank     = int(comm.rank());
    int mpi_size     = int(comm.size());

    auto xy = array::make_view<double, 2>(dynamic_cast<const functionspace::detail::StructuredColumns&>(fs_).xy());

    const std::vector<Point2>& points = config.getBool("inner_bounding_box", false) ? inner_bounding_box_ : points_;

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < points.size(); ++i) {
        xmin = std::min(xmin, points[i][XX]);
        xmax = std::max(xmax, points[i][XX]);
    }
    comm.allReduceInPlace(xmin, eckit::mpi::min());
    comm.allReduceInPlace(xmax, eckit::mpi::max());

    idx_t count     = dynamic_cast<const functionspace::detail::StructuredColumns&>(fs_).sizeOwned();
    idx_t count_all = count;
    comm.allReduceInPlace(count_all, eckit::mpi::sum());

    bool plot_nodes = config.getBool("nodes", false);
    for (int r = 0; r < mpi_size; ++r) {
        // clang-format off
        if ( mpi_rank == r ) {
            std::ofstream f( filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app );
            ATLAS_ASSERT( f.is_open() );
            if ( mpi_rank == 0 ) {
                f << "\n" "import sys"
                     "\n"
                     "\n" "# Configuration option to plot nodes"
                     "\n" "plot_nodes = False"
                     "\n" "for argv in sys.argv[1:] :"
                     "\n" "  if argv == \"--nodes\" :"
                     "\n" "      plot_nodes = " + std::string( plot_nodes ? "True" : "False" ) +
                     "\n"
                     "\n" "import matplotlib.pyplot as plt"
                     "\n" "from matplotlib.path import Path"
                     "\n" "import matplotlib.patches as patches"
                     "\n"
                     "\n" "from itertools import cycle"
                     "\n" "import matplotlib.cm as cm"
                     "\n" "import numpy as np"
                     "\n" "colours = cycle([cm.Paired(i) for i in np.linspace(0,1,12,endpoint=True)])"
                     "\n"
                     "\n" "fig = plt.figure()"
                     "\n" "ax = fig.add_subplot(111,aspect='equal')"
                     "\n";
            }
            if (points.size() > 0) {
                f << "\n" "verts_" << r << " = [";
                for ( size_t i=0; i<points.size(); ++i ) {
                    f << "\n  (" << points[i][XX] << ", " << points[i][YY] << "), ";
                }
                f << "\n]"
                    "\n"
                    "\n" "codes_" << r << " = [Path.MOVETO]"
                    "\n" "codes_" << r << ".extend([Path.LINETO] * " << ( points.size() - 2 ) << ")"
                    "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
                    "\n"
                    "\n" "count_" << r << " = " << count <<
                    "\n" "count_all_" << r << " = " << count_all <<
                    "\n";
                if ( plot_nodes ) {
                    f << "\n" "x_" << r << " = [";
                    for ( idx_t i = 0; i < count; ++i ) {
                        f << xy( i, XX ) << ", ";
                    }
                    f << "]"
                        "\n" "y_" << r << " = [";
                    for ( idx_t i = 0; i < count; ++i ) {
                        f << xy( i, YY ) << ", ";
                    }
                    f << "]";
                }
                f << "\n"
                    "\n" "c = next(colours)"
                    "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), edgecolor=c, facecolor=c, alpha=0.3, lw=1))";
                if ( plot_nodes ) {
                    f << "\n" "if plot_nodes:"
                        "\n" "    ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')";
                }
                f << "\n";
            }
            if ( mpi_rank == mpi_size - 1 ) {
                if (fs_.projection().units() == "degrees") {
                    f << "\n" "ax.set_xlim( " << xmin << "-5, " << xmax << "+5)"
                         "\n" "ax.set_ylim(-90-5,  90+5)"
                         "\n" "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
                         "\n" "ax.set_yticks([-90,-45,0,45,90])";
                }
                else {
                    f << "\n" "ax.autoscale()";
                }
                f << "\n" "plt.grid()"
                     "\n" "plt.show()";
            }
        }
        // clang-format on
        comm.barrier();
    }
}

util::PartitionPolygon::PointsXY StructuredPartitionPolygon::xy() const {
    return points_;
}

util::PartitionPolygon::PointsXY StructuredPartitionPolygon::lonlat() const {
    return points_;
}

void StructuredPartitionPolygon::allGather(util::PartitionPolygons& polygons_) const {
    ATLAS_TRACE();

    const auto& comm   = mpi::comm(fs_.mpi_comm());
    const int mpi_size = int(comm.size());

    polygons_.clear();
    polygons_.reserve(mpi_size);

    auto& poly = *this;

    std::vector<double> mypolygon;
    mypolygon.reserve(poly.size() * 2);

    for (auto& p : poly.xy()) {
        mypolygon.push_back(p[XX]);
        mypolygon.push_back(p[YY]);
    }

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
        polygons_.emplace_back(new util::ExplicitPartitionPolygon(std::move(recv_points)));
    }
}

// void StructuredPartitionPolygon::print( std::ostream& out ) const {
//     out << "polygon:{"
//         << "halo:" << halo_ << ",size:" << polygon_.size() << ",nodes:" << static_cast<const util::Polygon&>( polygon_ )
//         << "}";
// }

// ----------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
