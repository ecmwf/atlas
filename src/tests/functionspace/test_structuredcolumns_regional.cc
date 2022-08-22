/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"


#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


#include "tests/AtlasTestEnvironment.h"

using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("regional lonlat") {
    const int Nx = 8, Ny = 8;
    const double xmin = +20, xmax = +60, ymin = +20, ymax = +60;

    StructuredGrid::XSpace xspace(Config("type", "linear")("N", Nx)("start", xmin)("end", xmax));
    StructuredGrid::YSpace yspace(Config("type", "linear")("N", Ny)("start", ymin)("end", ymax));

    StructuredGrid grid(xspace, yspace);

    functionspace::StructuredColumns fs(grid, grid::Partitioner("checkerboard"), Config("halo", 1));

    ATLAS_TRACE_SCOPE("Output python") {
        auto xy = array::make_view<double, 2>(fs.xy());
        auto g  = array::make_view<gidx_t, 1>(fs.global_index());
        auto p  = array::make_view<int, 1>(fs.partition());


        eckit::PathName filepath("test_structuredcolumns_regional_p" + std::to_string(mpi::comm().rank()) + ".py");

        std::ofstream f(filepath.asString().c_str(), std::ios::trunc);

        f << "\n"
             "import matplotlib.pyplot as plt"
             "\n"
             "from matplotlib.path import Path"
             "\n"
             "import matplotlib.patches as patches"
             "\n"
             ""
             "\n"
             "from itertools import cycle"
             "\n"
             "import matplotlib.cm as cm"
             "\n"
             "import numpy as np"
             "\n"
             ""
             "\n"
             "fig = plt.figure(figsize=(20,10))"
             "\n"
             "ax = fig.add_subplot(111,aspect='equal')"
             "\n"
             "";

        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        f << "\n"
             "x = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << xy(n, XX) << ", ";
                xmin = std::min(xmin, xy(n, XX));
                xmax = std::max(xmax, xy(n, XX));
            }
        }
        f << "]";

        f << "\n"
             "y = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << xy(n, YY) << ", ";
                ymin = std::min(ymin, xy(n, YY));
                ymax = std::max(ymax, xy(n, YY));
            }
        }
        f << "]";

        f << "\n"
             "g = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << g(n) << ", ";
            }
        }
        f << "]";

        f << "\n"
             "p = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << p(n) << ", ";
            }
        }
        f << "]";

        //        f << "\n"
        //             "r = [";
        //        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        //            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
        //                idx_t n = fs.index(i, j);
        //                f << r(n) << ", ";
        //            }
        //        }
        //        f << "]";

        f << "\n"
             ""
             "\n"
             "c = [ cm.Paired( float(pp%13)/12. ) for pp in p ]"
             "\n"
             "ax.scatter(x, y, color=c, marker='o')"
             "\n"
             "for i in range("
          << fs.size()
          << "):"
             "\n"
             "  ax.annotate(g[i], (x[i],y[i]), fontsize=8)"
             "\n"
             "";
        f << "\n"
             "ax.set_xlim( "
          << std::min(0., xmin) << "-5, " << std::max(360., xmax)
          << "+5)"
             "\n"
             "ax.set_ylim( "
          << std::min(-90., ymin) << "-5, " << std::max(90., ymax)
          << "+5)"
             "\n"
             "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
             "\n"
             "ax.set_yticks([-90,-45,0,45,90])"
             "\n"
             "plt.grid()"
             "\n"
             "plt.show()"
             "\n";
    }

    auto lonlat = array::make_view<double, 2>(fs.lonlat());
    for (idx_t n = 0; n < lonlat.shape(0); ++n) {
        // Log::info() << n << "\t" << PointLonLat(lonlat(n, LON), lonlat(n, LAT)) << std::endl;
        EXPECT(lonlat(n, LON) >= xmin);
        EXPECT(lonlat(n, LON) <= xmax);
        EXPECT(lonlat(n, LAT) >= ymin);
        EXPECT(lonlat(n, LAT) <= ymax);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
