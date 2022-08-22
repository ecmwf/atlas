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
#include <cmath>
#include <iomanip>
#include <numeric>

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/local/TransLocal.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "tests/AtlasTestEnvironment.h"

#if ATLAS_HAVE_TRANS
#include "atlas/library/config.h"
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif
#endif

using namespace eckit;

using atlas::array::Array;
using atlas::array::ArrayView;
using atlas::array::make_view;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

struct AtlasTransEnvironment : public AtlasTestEnvironment {
    AtlasTransEnvironment(int argc, char* argv[]): AtlasTestEnvironment(argc, argv) {
#if ATLAS_HAVE_TRANS
        trans_use_mpi(mpi::comm().size() > 1);
        trans_init();
#endif
    }

#if ATLAS_HAVE_TRANS
    ~AtlasTransEnvironment() { trans_finalize(); }
#endif
};

//-----------------------------------------------------------------------------
// Routine to compute the spherical harmonics analytically at one point
// (up to wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
double sphericalharmonics_analytic_point(
    const int n,         // total wave number (implemented so far for n<4
    const int m,         // zonal wave number (implemented so far for m<4, m<n
    const int imag,      // 0: test real part, 1: test imaginary part
    const double lon,    // longitude in radians
    const double lat,    // latitude in radians
    const int ivar_in,   // variable that is set to 1 for wave number n,m. 0:
                         // vorticity, 1: divergence, 2: scalar
    const int ivar_out)  // variable returned by this function. 0: u, 1: v, 2: scalar
{
    double latsin = std::sin(lat), latcos = std::cos(lat);
    double lonsin = std::sin(m * lon), loncos = std::cos(m * lon);
    double a = util::Earth::radius();
    // Fourier part of the spherical harmonics:
    double rft = 1.;
    if (m > 0) {
        rft *= 2.;  // the famous factor 2 that noone really understands
    }
    if (imag == 0) {
        rft *= loncos;
    }
    else {
        rft *= -lonsin;
    }
    // Legendre part of the spherical harmonics (following
    // http://mathworld.wolfram.com/SphericalHarmonic.html
    // multiplied with -2*sqrt(pi) due to different normalization and different
    // coordinates):
    // (can also be computed on http://www.wolframalpha.com with:
    // LegendreP[n, m, x]/Sqrt[1/2*Integrate[LegendreP[n, m, y]^2, {y, -1, 1}]])
    // n, m need to be replaced by hand with the correct values
    // (otherwise the command will be too long for the free version of
    // wolframalpha)

    // scalar:
    if (ivar_in == 2) {
        if (ivar_out == 2) {
            if (m == 0 && n == 0) {
                return rft;
            }
            if (m == 0 && n == 1) {
                return std::sqrt(3.) * latsin * rft;
            }
            if (m == 0 && n == 2) {
                return std::sqrt(5.) / 2. * (3. * latsin * latsin - 1.) * rft;  // sign?
            }
            if (m == 0 && n == 3) {
                return std::sqrt(7.) / 2. * (5. * latsin * latsin - 3.) * latsin * rft;  // sign?
            }
            if (m == 1 && n == 1) {
                return std::sqrt(3. / 2.) * latcos * rft;  // sign?
            }
            if (m == 1 && n == 2) {
                return std::sqrt(15. / 2.) * latsin * latcos * rft;  // sign?
            }
            if (m == 1 && n == 3) {
                return std::sqrt(21.) / 4. * latcos * (5. * latsin * latsin - 1.) * rft;  // sign?
            }
            if (m == 2 && n == 2) {
                return std::sqrt(15. / 2.) / 2. * latcos * latcos * rft;
            }
            if (m == 2 && n == 3) {
                return std::sqrt(105. / 2.) / 2. * latcos * latcos * latsin * rft;
            }
            if (m == 3 && n == 3) {
                return std::sqrt(35.) / 4. * latcos * latcos * latcos * rft;  // sign?
            }
            if (m == 4 && n == 4) {
                return (3 * std::sqrt(17.5) * std::pow(latcos, 4)) / 8. * rft;
            }
            if (m == 5 && n == 5) {
                return (3 * std::sqrt(77) * std::pow(latcos, 5)) / 16. * rft;
            }
            if (m == 6 && n == 6) {
                return (std::sqrt(3003) * std::pow(latcos, 6)) / 32. * rft;
            }
            if (m == 7 && n == 7) {
                return (3 * std::sqrt(357.5) * std::pow(latcos, 7)) / 32. * rft;
            }
            if (m == 8 && n == 8) {
                return (3 * std::sqrt(6077.5) * std::pow(latcos, 8)) / 128. * rft;
            }
            if (m == 9 && n == 9) {
                return (std::sqrt(230945) * std::pow(latcos, 9)) / 256. * rft;
            }
            if (m == 10 && n == 10) {
                return (std::sqrt(969969) * std::pow(latcos, 10)) / 512. * rft;
            }
            if (m == 11 && n == 11) {
                return (std::sqrt(1.0140585e6) * std::pow(latcos, 11)) / 512. * rft;
            }
            if (m == 12 && n == 12) {
                return (5 * std::sqrt(676039) * std::pow(latcos, 12)) / 2048. * rft;
            }
            if (m == 13 && n == 13) {
                return (15 * std::sqrt(78004.5) * std::pow(latcos, 13)) / 2048. * rft;
            }
            if (m == 14 && n == 14) {
                return (15 * std::sqrt(323161.5) * std::pow(latcos, 14)) / 4096. * rft;
            }
            if (m == 15 && n == 15) {
                return (3 * std::sqrt(33393355) * std::pow(latcos, 15)) / 8192. * rft;
            }
            if (m == 16 && n == 16) {
                return (3 * std::sqrt(5.509903575e8) * std::pow(latcos, 16)) / 32768. * rft;
            }
            if (m == 17 && n == 17) {
                return (15 * std::sqrt(90751353) * std::pow(latcos, 17)) / 65536. * rft;
            }
            if (m == 18 && n == 18) {
                return (5 * std::sqrt(3357800061) * std::pow(latcos, 18)) / 131072. * rft;
            }
            if (m == 19 && n == 19) {
                return (15 * std::sqrt(3.829070245e8) * std::pow(latcos, 19)) / 131072. * rft;
            }
            if (m == 20 && n == 20) {
                return (3 * std::sqrt(156991880045) * std::pow(latcos, 20)) / 524288. * rft;
            }
            if (m == 21 && n == 21) {
                return (std::sqrt(1.4465680375575e12) * std::pow(latcos, 21)) / 524288. * rft;
            }
            if (m == 22 && n == 22) {
                return (15 * std::sqrt(2.63012370465e10) * std::pow(latcos, 22)) / 1.048576e6 * rft;
            }
            if (m == 23 && n == 23) {
                return (15 * std::sqrt(107492012277) * std::pow(latcos, 23)) / 2.097152e6 * rft;
            }
            if (m == 24 && n == 24) {
                return (105 * std::sqrt(35830670759) * std::pow(latcos, 24)) / 8.388608e6 * rft;
            }
            if (m == 25 && n == 25) {
                return (21 * std::sqrt(9.136821043545e11) * std::pow(latcos, 25)) / 8.388608e6 * rft;
            }
            if (m == 26 && n == 26) {
                return (21 * std::sqrt(3.7250116562145e12) * std::pow(latcos, 26)) / 1.6777216e7 * rft;
            }
            if (m == 27 && n == 27) {
                return (7 * std::sqrt(136583760727865.) * std::pow(latcos, 27)) / 3.3554432e7 * rft;
            }
            if (m == 28 && n == 28) {
                return (std::sqrt(2.7248460265209068e16) * std::pow(latcos, 28)) / 6.7108864e7 * rft;
            }
            if (m == 29 && n == 29) {
                return (std::sqrt(110873045217057585.) * std::pow(latcos, 29)) / 1.34217728e8 * rft;
            }
            if (m == 30 && n == 30) {
                return (std::sqrt(450883717216034179.) * std::pow(latcos, 30)) / 2.68435456e8 * rft;
            }
            if (m == 31 && n == 31) {
                return (21 * std::sqrt(1.0389025742304935e15) * std::pow(latcos, 31)) / 2.68435456e8 * rft;
            }
            if (m == 32 && n == 32) {
                return (21 * std::sqrt(6.752866732498208e16) * std::pow(latcos, 32)) / 2.147483648e9 * rft;
            }
            if (m == 33 && n == 33) {
                return (7 * std::sqrt(2467865842240254105.) * std::pow(latcos, 33)) / 4.294967296e9 * rft;
            }
            if (m == 34 && n == 34) {
                return (21 * std::sqrt(1112959105324036165.) * std::pow(latcos, 34)) / 8.589934592e9 * rft;
            }
            if (m == 35 && n == 35) {
                return (3 * std::sqrt(5.53140675346046e19) * std::pow(latcos, 35)) / 8.589934592e9 * rft;
            }
            if (m == 36 && n == 36) {
                return (std::sqrt(8075853860052271220473.) * std::pow(latcos, 36)) / 3.4359738368e10 * rft;
            }
            if (m == 37 && n == 37) {
                return (5 * std::sqrt(3.2739948081292994e20) * std::pow(latcos, 37)) / 3.4359738368e10 * rft;
            }
            if (m == 38 && n == 38) {
                return (35 * std::sqrt(2.707815254843781e19) * std::pow(latcos, 38)) / 6.8719476736e10 * rft;
            }
            if (m == 39 && n == 39) {
                return (35 * std::sqrt(109701233401363445369.) * std::pow(latcos, 39)) / 1.37438953472e11 * rft;
            }
            if (m == 40 && n == 40) {
                return (63 * std::sqrt(548506167006817226845.) * std::pow(latcos, 40)) / 5.49755813888e11 * rft;
            }
            if (m == 41 && n == 41) {
                return (63 * std::sqrt(5.551952666044613e20) * std::pow(latcos, 41)) / 5.49755813888e11 * rft;
            }
            if (m == 42 && n == 42) {
                return (15 * std::sqrt(3.964094203555854e22) * std::pow(latcos, 42)) / 1.099511627776e12 * rft;
            }
            if (m == 43 && n == 43) {
                return (45 * std::sqrt(17823059209786010066579.) * std::pow(latcos, 43)) / 2.199023255552e12 * rft;
            }
            if (m == 44 && n == 44) {
                return (45 * std::sqrt(7.210237589413431e22) * std::pow(latcos, 44)) / 4.398046511104e12 * rft;
            }
            if (m == 45 && n == 45) {
                return (21 * std::sqrt(1339044123748208678378695.) * std::pow(latcos, 45)) / 8.796093022208e12 * rft;
            }
            //    return
            //    std::pow(latcos,45)*rft*21.*std::sqrt(1339044123748208678378695.)/8796093022208.;
            //    // sign?
        }
        else {
            return 0.;
        }
    }

    // for the remainder the factor 2 from rft is already included in the
    // formulas:

    // vorticity:
    if (ivar_in == 0) {
        if (ivar_out == 0) {  // u:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                if (imag == 0) {
                    return std::sqrt(3.) * a / 2. * latcos;
                }
                else {
                    return 0.;
                }
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return -a * std::sqrt(3. / 2.) * loncos * latsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * lonsin * latsin;
                }
            }
        }
        else if (ivar_out == 1) {  // v:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                return 0.;
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * lonsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * loncos;
                }
            }
        }
        else {
            return 0.;
        }
    }

    // divergence:
    if (ivar_in == 1) {
        if (ivar_out == 0) {  // u:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                return 0.;
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * lonsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * loncos;
                }
            }
        }
        else if (ivar_out == 1) {  // v:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                if (imag == 0) {
                    return -std::sqrt(3.) * a / 2. * latcos;
                }
                else {
                    return 0.;
                }
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * loncos * latsin;
                }
                else {
                    return -a * std::sqrt(3. / 2.) * lonsin * latsin;  // sign?
                }
            }
        }
        else {
            return 0.;
        }
    }

    return -1.;
}

//-----------------------------------------------------------------------------
// Routine to compute the spherical harmonics analytically on a grid
// (up to wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
void spectral_transform_grid_analytic(const int trc,       // truncation (in)
                                      bool trcFT,          // truncation for Fourier transformation (in)
                                      const int n,         // total wave number (implemented so far for n<4)
                                      const int m,         // zonal wave number (implemented so far for m<4, m<n)
                                      const int imag,      // 0: test real part, 1: test imaginary part
                                      const Grid grid,     // call with something like Grid("O32")
                                      double rspecg[],     // spectral data, size (trc+1)*trc (out)
                                      double rgp[],        // resulting grid point data (out)
                                      const int ivar_in,   // variable that is set to 1 for wave number n,m. 0:
                                                           // vorticity, 1: divergence, 2: scalar
                                      const int ivar_out)  // variable returned by this function. 0: u, 1: v, 2: scalar
{
    int N = (trc + 2) * (trc + 1) / 2;
    for (int jm = 0; jm < 2 * N; jm++) {
        rspecg[jm] = 0.;
    }
    int k = 0;
    for (int jm = 0; jm <= trc; jm++) {
        for (int jn = jm; jn <= trc; jn++) {
            if (jm == m && jn == n) {
                rspecg[2 * k + imag]       = 1.;
                rspecg[2 * k + (1 - imag)] = 0.;
            }
            k++;
        }
    }

    for (int jm = 0; jm < grid.size(); jm++) {
        rgp[jm] = 0.;
    }

    if (StructuredGrid(grid)) {
        StructuredGrid g(grid);
        Grid gridGlobal;
        StructuredGrid gs_global;
        int jlatMin = 0;
        if (trcFT) {
            gridGlobal      = Grid(grid.name());
            gs_global       = StructuredGrid(gridGlobal);
            int nlatsGlobal = gs_global.ny();
            for (int jlat = 0; jlat < nlatsGlobal; jlat++) {
                if (gs_global.y(jlat) > g.y(0)) {
                    jlatMin++;
                }
            }
        }

        int idx = 0;
        for (idx_t j = 0; j < g.ny(); ++j) {
            double lat = g.y(j) * util::Constants::degreesToRadians();
            int ftrc   = trc + 1;
            if (trcFT) {
                ftrc = trans::fourier_truncation(trc, gs_global.nx(jlatMin + j), gs_global.nxmax(), gs_global.ny(), lat,
                                                 RegularGrid(gs_global));
            }
            /*Log::info() << "j=" << j << " ftrc=" << ftrc << " trc=" << trc << " nx=" << gs_global.nx( jlatMin + j )
                        << " nxmax=" << gs_global.nxmax() << " nlats=" << gs_global.ny() << " lat=" << g.y( j )
                        << " jlatMin=" << jlatMin << std::endl;*/
            for (idx_t i = 0; i < g.nx(j); ++i) {
                double lon = g.x(i, j) * util::Constants::degreesToRadians();

                // compute spherical harmonics:
                if (ftrc > m) {
                    rgp[idx++] = sphericalharmonics_analytic_point(n, m, imag, lon, lat, ivar_in, ivar_out);
                }
                else {
                    rgp[idx++] = 0.;
                }
            }
        }
    }
    else {
        int idx = 0;
        for (const PointXY& p : grid.xy()) {
            double lon = p.x() * util::Constants::degreesToRadians();
            double lat = p.y() * util::Constants::degreesToRadians();
            // compute spherical harmonics:
            rgp[idx++] = sphericalharmonics_analytic_point(n, m, imag, lon, lat, ivar_in, ivar_out);
        }
    }
}

//-----------------------------------------------------------------------------
// Compute root mean square of difference between two arrays divided by maximum
// absolute value of second array
//
// Author:
// Andreas Mueller *ECMWF*
//
double compute_rms(const size_t N,   // length of the arrays
                   double array1[],  // first of the two arrays
                   double array2[])  // second of the two arrays
{
    double rms = 0., rmax = 0.;
    for (size_t idx = 0; idx < N; idx++) {
        double diff = array1[idx] - array2[idx];
        rms += diff * diff;
        rmax = std::max(rmax, std::abs(array2[idx]));
    }
    if (rmax == 0.) {
        rms = 0.;
    }
    else {
        rms = std::sqrt(rms / N) / rmax;
    }
    return rms;
}

//-----------------------------------------------------------------------------
#if 1
CASE("test_trans_vordiv_with_translib") {
    Log::info() << "test_trans_vordiv_with_translib" << std::endl;
    // test transgeneral by comparing its result with the trans library
    // this test is based on the test_nomesh case in test_trans.cc

    double tolerance = 1.e-13;

    // Grid: (Adjust the following line if the test takes too long!)
    Grid g("F64");

    StructuredGrid gs(g);
    int ndgl = gs.ny();
    //int trc  = ndgl - 1;  // linear
    int trc = ndgl / 2. - 1;  // cubic
#if ATLAS_HAVE_TRANS
    trans::Trans transIFS(g, trc, util::Config("type", "ectrans"));
    double rav = 0.;  // compute average rms error of trans library in rav
#endif
    trans::Trans transLocal1(g, trc, util::Config("type", "local"));
    trans::Trans transLocal2(g, trc, util::Config("type", "local"));
    double rav1 = 0., rav2 = 0.;  // compute average rms errors of transLocal1 and transLocal2

    functionspace::Spectral spectral(trc);
    functionspace::StructuredColumns gridpoints(g);

    int nb_scalar = 1, nb_vordiv = 0;
    int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
    std::vector<double> sp(2 * N * nb_scalar);
    std::vector<double> vor(2 * N * nb_vordiv);
    std::vector<double> div(2 * N * nb_vordiv);
    std::vector<double> rspecg(2 * N);
    std::vector<double> gp(nb_all * g.size());
    std::vector<double> rgp1(nb_all * g.size());
    std::vector<double> rgp2(nb_all * g.size());
    std::vector<double> rgp_analytic(g.size());

    int icase = 0;
    for (int ivar_in = 2; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
        for (int ivar_out = 2; ivar_out < 3; ivar_out++) {  // u, v, scalar
            int nb_fld = 1;
            if (ivar_out == 2) {
                tolerance = 1.e-13;
                nb_fld    = nb_scalar;
            }
            else {
                tolerance = 2.e-6;
                nb_fld    = nb_vordiv;
            }
            for (int jfld = 0; jfld < nb_fld; jfld++) {  // multiple fields
                int k = 0;
                for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                    for (int n = m; n <= trc; n++) {             // total wavenumber
                        for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                            if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0.) {
                                for (int j = 0; j < 2 * N * nb_scalar; j++) {
                                    sp[j] = 0.;
                                }
                                for (int j = 0; j < 2 * N * nb_vordiv; j++) {
                                    vor[j] = 0.;
                                    div[j] = 0.;
                                }
                                if (ivar_in == 0) {
                                    vor[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 1) {
                                    div[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 2) {
                                    sp[k * nb_scalar + jfld] = 1.;
                                }

                                for (int j = 0; j < nb_all * g.size(); j++) {
                                    gp[j]   = 0.;
                                    rgp1[j] = 0.;
                                    rgp2[j] = 0.;
                                }
                                for (int j = 0; j < g.size(); j++) {
                                    rgp_analytic[j] = 0.;
                                }

                                spectral_transform_grid_analytic(trc, trc, n, m, imag, g, rspecg.data(),
                                                                 rgp_analytic.data(), ivar_in, ivar_out);

                                EXPECT_NO_THROW(transLocal1.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp1.data()));

                                EXPECT_NO_THROW(transLocal2.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp2.data()));

                                int pos = (ivar_out * nb_vordiv + jfld);

                                double rms_gen1 =
                                    compute_rms(g.size(), rgp1.data() + pos * g.size(), rgp_analytic.data());

                                double rms_gen2 =
                                    compute_rms(g.size(), rgp2.data() + pos * g.size(), rgp_analytic.data());

                                rav1 += rms_gen1;
                                rav2 += rms_gen2;
                                if (!(rms_gen1 < tolerance) || !(rms_gen2 < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
                                EXPECT(rms_gen1 < tolerance);
                                EXPECT(rms_gen2 < tolerance);
                                icase++;

#if ATLAS_HAVE_TRANS
                                EXPECT_NO_THROW(transIFS.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                  div.data(), gp.data()));
                                double rms_trans =
                                    compute_rms(g.size(), gp.data() + pos * g.size(), rgp_analytic.data());
                                rav += rms_trans;
                                double rms_diff =
                                    compute_rms(g.size(), rgp1.data() + pos * g.size(), gp.data() + pos * g.size());
                                EXPECT(rms_trans < tolerance);
                                if (!(rms_trans < tolerance) || !(rms_diff < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(rms_trans);
                                    ATLAS_DEBUG_VAR(rms_diff);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
#endif
                                EXPECT(icase < 300);
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }
    Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
    rav1 /= icase;
    Log::info() << "average RMS error of transLocal1: " << rav1 << std::endl;
    rav2 /= icase;
    Log::info() << "average RMS error of transLocal2: " << rav2 << std::endl;
#if ATLAS_HAVE_TRANS
    rav /= icase;
    Log::info() << "average RMS error of transIFS: " << rav << std::endl;
#endif
}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE("test_trans_hires") {
    Log::info() << "test_trans_hires" << std::endl;
    // test transgeneral by comparing its result with the trans library
    // this test is based on the test_nomesh case in test_trans.cc

    std::ostream& out = Log::info();
    double tolerance  = 1.e-13;

    //#if ATLAS_HAVE_TRANS
    //    //std::string transTypes[4] = {"localopt", "localopt2", "Local", "ectrans"};
    //    //std::string transTypes[2] = {"localopt2", "Local"};
    //    //std::string transTypes[3] = {"Local", "localopt2", "localopt"};
    //    //std::string transTypes[1] = {"Local"};
    //#else
    //    //std::string transTypes[1] = {"localopt2"};
    //#endif
    std::vector<util::Config> trans_configs;
    std::vector<std::string> transTypes = {"local"};
    if (trans::Trans::hasBackend("ectrans")) {
        transTypes.emplace_back("ectrans");
    }
    std::map<std::string, std::vector<std::string>> backends{{"ectrans", {"lapack"}},
                                                             {"local", {"generic", "openmp", "eigen"}}};

    //Domain testdomain =F ZonalBandDomain( {-90., 90.} );
    //Domain testdomain = ZonalBandDomain( {-.5, .5} );
    //Domain testdomain = RectangularDomain( {0., 30.}, {-.05, .05} );
    //Domain testdomain = ZonalBandDomain( {-85., -86.} );
    //Domain testdomain = RectangularDomain( {-.01, .01}, {-.01, .01} );
    //Domain testdomain = RectangularDomain( {-1, 1}, {-1, 1} );
    Domain testdomain = GlobalDomain();
    // Grid: (Adjust the following line if the test takes too long!)
    Grid g("O160", testdomain);
    Grid g_global(g.name());

    StructuredGrid gs(g);
    StructuredGrid gs_global(g_global);
    Log::info() << "nlats: " << gs.ny() << " nlons:" << gs.nxmax() << std::endl;
    int ndgl = gs_global.ny();
    //int trc  = ndgl - 1;  // linear
    int trc = ndgl / 2. - 1;  // cubic
    Log::info() << "truncation: " << trc << std::endl;

    int nb_scalar = 1, nb_vordiv = 0;

    for (auto& trans_type : transTypes) {
        ATLAS_TRACE(trans_type);
        int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
        int icase = 0;
        trans::Trans trans(g, trc, option::type(trans_type));
        for (int ivar_in = 2; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
            for (int ivar_out = 2; ivar_out < 3; ivar_out++) {  // u, v, scalar
                int nb_fld = 1;
                if (ivar_out == 2) {
                    nb_fld = nb_scalar;
                }
                else {
                    nb_fld = nb_vordiv;
                }
                for (int jfld = 0; jfld < 1; jfld++) {  // multiple fields
                    int k = 0;
                    for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                        for (int n = m; n <= trc; n++) {             // total wavenumber
                            for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                                if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. &&
                                    icase < 1) {
                                    std::vector<double> sp(2 * N * nb_scalar);
                                    std::vector<double> gp(nb_all * g.size());
                                    if (ivar_in == 2) {
                                        sp[k * nb_scalar + jfld] = 1.;
                                    }
                                    auto original_backend = linalg::dense::current_backend();
                                    for (auto& backend : backends[trans_type]) {
                                        if (linalg::dense::Backend(backend).available()) {
                                            linalg::dense::current_backend(backend);
                                            ATLAS_TRACE("invtrans [backend=" + backend + "]");
                                            auto start = std::chrono::system_clock::now();
                                            EXPECT_NO_THROW(trans.invtrans(nb_scalar, sp.data(), nb_vordiv, nullptr,
                                                                           nullptr, gp.data()));
                                            auto end = std::chrono::system_clock::now();  //
                                            std::chrono::duration<double> elapsed_seconds = end - start;
                                            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                                            std::string time_str = std::ctime(&end_time);
                                            Log::info()
                                                << trans_type << "[backend=" << backend << "] : case " << icase
                                                << ", elapsed time: " << elapsed_seconds.count()
                                                << "s. Now: " << time_str.substr(0, time_str.length() - 1) << std::endl;
                                        }
                                    }
                                    linalg::dense::current_backend(original_backend);
                                    icase++;
                                }
                                k++;
                            }
                        }
                    }
                }
            }
        }
        //.Log::info() << "Vordiv+scalar comparison with trans::" << description << ": all " << icase
        //            << " cases successfully passed!" << std::endl;
    }
}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE("test_trans_domain") {
    Log::info() << "test_trans_domain" << std::endl;
    // test transgeneral by comparing with analytic solution on a cropped domain
    // this test also includes testing caching

    //Domain testdomain = ZonalBandDomain( {-90., 90.} );
    //Domain testdomain = ZonalBandDomain( {-.5, .5} );
    //Domain testdomain = RectangularDomain( {0., 30.}, {-.05, .05} );
    //Domain testdomain1 = ZonalBandDomain( {-10., 5.} );
    Domain testdomain1 = RectangularDomain({-5., 5.}, {-2.5, 0.});
    //Domain testdomain1 = RectangularDomain( {-1., 1.}, {50., 55.} );
    Domain testdomain2 = RectangularDomain({-5., 5.}, {-2.5, 0.});
    // Grid: (Adjust the following line if the test takes too long!)

    Grid global_grid("O64");
    //Grid g1( global_grid, testdomain1 );
    Grid g2(global_grid, testdomain2);
    Grid g1(global_grid);
    //Grid g2( global_grid );

    bool fourierTrc1 = true;
    bool fourierTrc2 = false;
    // fourierTrc1, fourierTrc2: need to be false if no global grid can be constructed
    //                           (like for grids created with LinearSpacing)

    //using LinearSpacing = grid::LinearSpacing;
    //StructuredGrid g2( LinearSpacing( {0., 180.}, 3 ), LinearSpacing( {89., 90.}, 2 ) );
    // when using LinearSpacing: set fourierTrc2 to false

    int trc = 63;
    //Log::info() << "rgp1:" << std::endl;
    if (eckit::PathName("legcache.bin").exists()) {
        eckit::PathName("legcache.bin").unlink();
    }
    Trace t1(Here(), "translocal1 construction");
    trans::Trans transLocal1(global_grid, g1.domain(), trc,
                             option::type("local") | option::write_legendre("legcache.bin"));
    t1.stop();
    //Log::info() << "rgp2:" << std::endl;
    trans::Cache cache;
    ATLAS_TRACE_SCOPE("Read cache") cache = trans::LegendreCache("legcache.bin");
    Trace t2(Here(), "translocal2 construction");
    trans::Trans transLocal2(cache, global_grid, g2.domain(), trc,
                             option::type("local") | option::write_legendre("legcache2.bin"));
    //trans::Trans transLocal2( cache, g2, trc, option::type( "local" ) );
    //trans::Trans transLocal2( cache, g2, trc,
    //                          option::type( "local" ) | option::no_fft() );
    //trans::Trans transLocal2( g2, trc, option::type( "local" ) );
    t2.stop();

    double rav1 = 0., rav2 = 0.;  // compute average rms errors of transLocal1 and transLocal2

    functionspace::Spectral spectral(trc);

    int nb_scalar = 1, nb_vordiv = 1;
    int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
    std::vector<double> sp(2 * N * nb_scalar);
    std::vector<double> vor(2 * N * nb_vordiv);
    std::vector<double> div(2 * N * nb_vordiv);
    std::vector<double> rspecg(2 * N);
    std::vector<double> rgp1(nb_all * g1.size());
    std::vector<double> rgp2(nb_all * g2.size());
    std::vector<double> rgp1_analytic(g1.size());
    std::vector<double> rgp2_analytic(g2.size());

    StructuredGrid gs1(g1);
    StructuredGrid gs2(g2);
    double latPole   = 89.9999999;
    bool gridHasPole = gs1.y(0) > latPole || gs2.y(0) > latPole;
    double tolerance;
    int icase = 0;
    for (int ivar_in = 0; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
        for (int ivar_out = 0; ivar_out < 3; ivar_out++) {  // u, v, scalar
            int nb_fld = 1;
            if (ivar_out == 2 && !gridHasPole) {
                tolerance = 1.e-13;
                nb_fld    = nb_scalar;
            }
            else {
                tolerance = 2.e-6;
                nb_fld    = nb_vordiv;
            }
            for (int jfld = 0; jfld < nb_fld; jfld++) {  // multiple fields
                int k = 0;
                for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                    for (int n = m; n <= trc; n++) {             // total wavenumber
                        for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                            if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. &&
                                icase < 1000) {
                                auto start = std::chrono::system_clock::now();
                                for (int j = 0; j < 2 * N * nb_scalar; j++) {
                                    sp[j] = 0.;
                                }
                                for (int j = 0; j < 2 * N * nb_vordiv; j++) {
                                    vor[j] = 0.;
                                    div[j] = 0.;
                                }
                                if (ivar_in == 0) {
                                    vor[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 1) {
                                    div[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 2) {
                                    sp[k * nb_scalar + jfld] = 1.;
                                }

                                for (int j = 0; j < nb_all * g1.size(); j++) {
                                    rgp1[j] = 0.;
                                }
                                for (int j = 0; j < nb_all * g2.size(); j++) {
                                    rgp2[j] = 0.;
                                }
                                for (int j = 0; j < g1.size(); j++) {
                                    rgp1_analytic[j] = 0.;
                                }
                                for (int j = 0; j < g2.size(); j++) {
                                    rgp2_analytic[j] = 0.;
                                }

                                spectral_transform_grid_analytic(trc, fourierTrc1, n, m, imag, g1, rspecg.data(),
                                                                 rgp1_analytic.data(), ivar_in, ivar_out);

                                spectral_transform_grid_analytic(trc, fourierTrc2, n, m, imag, g2, rspecg.data(),
                                                                 rgp2_analytic.data(), ivar_in, ivar_out);

                                //Log::info() << std::endl << "rgp1:";
                                ATLAS_TRACE_SCOPE("translocal1")
                                EXPECT_NO_THROW(transLocal1.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp1.data()));

                                //Log::info() << std::endl << "rgp2:";
                                ATLAS_TRACE_SCOPE("translocal2")
                                EXPECT_NO_THROW(transLocal2.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp2.data()));

                                int pos = (ivar_out * nb_vordiv + jfld);

                                double rms_gen1 =
                                    compute_rms(g1.size(), rgp1.data() + pos * g1.size(), rgp1_analytic.data());

                                double rms_gen2 =
                                    compute_rms(g2.size(), rgp2.data() + pos * g2.size(), rgp2_analytic.data());

                                //Log::info() << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                //            << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl
                                //            << "rgp1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1[pos * g1.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl << "rgp2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2[pos * g2.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl;
                                rav1 += rms_gen1;
                                rav2 += rms_gen2;
                                if (!(rms_gen1 < tolerance) || !(rms_gen2 < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
                                EXPECT(rms_gen1 < tolerance);
                                EXPECT(rms_gen2 < tolerance);
                                icase++;
                                auto end                                      = std::chrono::system_clock::now();  //
                                std::chrono::duration<double> elapsed_seconds = end - start;
                                std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                                std::string time_str = std::ctime(&end_time);
                                Log::debug() << "case " << icase << ", elapsed time: " << elapsed_seconds.count()
                                             << "s. Now: " << time_str.substr(0, time_str.length() - 1) << std::endl;
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }
    Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
    rav1 /= icase;
    Log::info() << "average RMS error of transLocal1: " << rav1 << std::endl;
    rav2 /= icase;
    Log::info() << "average RMS error of transLocal2: " << rav2 << std::endl;
}
//-----------------------------------------------------------------------------
CASE("test_trans_pole") {
    Log::info() << "test_trans_pole" << std::endl;
    // test transform at the pole and with LinearSpacing grids
    // not using caching in this test because not useful for LinearSpacing grids

    //Domain testdomain = ZonalBandDomain( {-90., 90.} );
    //Domain testdomain = ZonalBandDomain( {-.5, .5} );
    //Domain testdomain = RectangularDomain( {0., 30.}, {-.05, .05} );
    //Domain testdomain1 = ZonalBandDomain( {-10., 5.} );
    Domain testdomain1 = RectangularDomain({-5., 5.}, {-2.5, 0.});
    //Domain testdomain1 = RectangularDomain( {-1., 1.}, {50., 55.} );
    Domain testdomain2 = RectangularDomain({-5., 5.}, {-2.5, 0.});
    // Grid: (Adjust the following line if the test takes too long!)

    Grid global_grid1("L3");
    Grid global_grid2("L3");
    //Grid g1( global_grid, testdomain1 );
    //Grid g2( global_grid, testdomain2 );
    Grid g1(global_grid1);
    //Grid g2( global_grid2 );

    bool fourierTrc1 = true;
    bool fourierTrc2 = false;
    // fourierTrc1, fourierTrc2: need to be false if no global grid can be constructed
    //                           (like for grids created with LinearSpacing)

    using LinearSpacing = grid::LinearSpacing;
    StructuredGrid g2(LinearSpacing({0., 180.}, 3), LinearSpacing({89., 90.}, 2));
    // when using LinearSpacing: set fourierTrc2 to false

    int trc = 2;
    Trace t1(Here(), "translocal1 construction");
    trans::Trans transLocal1(global_grid1, g1.domain(), trc, option::type("local"));
    t1.stop();
    Trace t2(Here(), "translocal2 construction");
    trans::Trans transLocal2(g2, trc, option::type("local"));
    t2.stop();

    double rav1 = 0., rav2 = 0.;  // compute average rms errors of transLocal1 and transLocal2

    functionspace::Spectral spectral(trc);

    int nb_scalar = 1, nb_vordiv = 1;
    int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
    std::vector<double> sp(2 * N * nb_scalar);
    std::vector<double> vor(2 * N * nb_vordiv);
    std::vector<double> div(2 * N * nb_vordiv);
    std::vector<double> rspecg(2 * N);
    std::vector<double> rgp1(nb_all * g1.size());
    std::vector<double> rgp2(nb_all * g2.size());
    std::vector<double> rgp1_analytic(g1.size());
    std::vector<double> rgp2_analytic(g2.size());

    StructuredGrid gs1(g1);
    StructuredGrid gs2(g2);
    double latPole   = 89.9999999;
    bool gridHasPole = gs1.y(0) > latPole || gs2.y(0) > latPole;
    double tolerance;
    int icase = 0;
    for (int ivar_in = 0; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
        for (int ivar_out = 0; ivar_out < 3; ivar_out++) {  // u, v, scalar
            int nb_fld = 1;
            if (ivar_out == 2 && !gridHasPole) {
                tolerance = 1.e-13;
                nb_fld    = nb_scalar;
            }
            else {
                tolerance = 2.e-5;
                nb_fld    = nb_vordiv;
            }
            for (int jfld = 0; jfld < nb_fld; jfld++) {  // multiple fields
                int k = 0;
                for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                    for (int n = m; n <= trc; n++) {             // total wavenumber
                        for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                            if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. &&
                                icase < 1000) {
                                auto start = std::chrono::system_clock::now();
                                for (int j = 0; j < 2 * N * nb_scalar; j++) {
                                    sp[j] = 0.;
                                }
                                for (int j = 0; j < 2 * N * nb_vordiv; j++) {
                                    vor[j] = 0.;
                                    div[j] = 0.;
                                }
                                if (ivar_in == 0) {
                                    vor[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 1) {
                                    div[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 2) {
                                    sp[k * nb_scalar + jfld] = 1.;
                                }

                                for (int j = 0; j < nb_all * g1.size(); j++) {
                                    rgp1[j] = 0.;
                                }
                                for (int j = 0; j < nb_all * g2.size(); j++) {
                                    rgp2[j] = 0.;
                                }
                                for (int j = 0; j < g1.size(); j++) {
                                    rgp1_analytic[j] = 0.;
                                }
                                for (int j = 0; j < g2.size(); j++) {
                                    rgp2_analytic[j] = 0.;
                                }

                                spectral_transform_grid_analytic(trc, fourierTrc1, n, m, imag, g1, rspecg.data(),
                                                                 rgp1_analytic.data(), ivar_in, ivar_out);

                                spectral_transform_grid_analytic(trc, fourierTrc2, n, m, imag, g2, rspecg.data(),
                                                                 rgp2_analytic.data(), ivar_in, ivar_out);

                                //Log::info() << std::endl << "rgp1:";
                                ATLAS_TRACE_SCOPE("translocal1")
                                EXPECT_NO_THROW(transLocal1.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp1.data()));

                                //Log::info() << std::endl << "rgp2:";
                                ATLAS_TRACE_SCOPE("translocal2")
                                EXPECT_NO_THROW(transLocal2.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp2.data()));

                                int pos = (ivar_out * nb_vordiv + jfld);

                                double rms_gen1 =
                                    compute_rms(g1.size(), rgp1.data() + pos * g1.size(), rgp1_analytic.data());

                                double rms_gen2 =
                                    compute_rms(g2.size(), rgp2.data() + pos * g2.size(), rgp2_analytic.data());

                                //Log::info() << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                //            << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl
                                //            << "rgp1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1[pos * g1.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl << "rgp2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2[pos * g2.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl;
                                rav1 += rms_gen1;
                                rav2 += rms_gen2;
                                if (!(rms_gen1 < tolerance) || !(rms_gen2 < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
                                EXPECT(rms_gen1 < tolerance);
                                EXPECT(rms_gen2 < tolerance);
                                icase++;
                                auto end                                      = std::chrono::system_clock::now();  //
                                std::chrono::duration<double> elapsed_seconds = end - start;
                                std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                                std::string time_str = std::ctime(&end_time);
                                Log::debug() << "case " << icase << ", elapsed time: " << elapsed_seconds.count()
                                             << "s. Now: " << time_str.substr(0, time_str.length() - 1) << std::endl;
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }
    Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
    rav1 /= icase;
    Log::info() << "average RMS error of transLocal1: " << rav1 << std::endl;
    rav2 /= icase;
    Log::info() << "average RMS error of transLocal2: " << rav2 << std::endl;
}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE("test_trans_southpole") {
    Log::info() << "test_trans_southpole" << std::endl;
    // test created for MIR-283 (limited area domain on the southern hemisphere with L-grid)


    //Domain testdomain = ZonalBandDomain( {-90., 90.} );
    //Domain testdomain = ZonalBandDomain( {-.5, .5} );
    //Domain testdomain = RectangularDomain( {0., 30.}, {-.05, .05} );
    //Domain testdomain1 = ZonalBandDomain( {-10., 5.} );
    //Domain testdomain1 = RectangularDomain( {-5., 5.}, {-2.5, 0.} );
    Domain testdomain1 = RectangularDomain({0., 10.}, {-90., -10.});
    //Domain testdomain1 = RectangularDomain( {-1., 1.}, {50., 55.} );
    //Domain testdomain2 = RectangularDomain( {-5., 5.}, {-2.5, 0.} );
    Domain testdomain2 = RectangularDomain({0., 10.}, {10., 90.});
    // Grid: (Adjust the following line if the test takes too long!)

    Grid global_grid1("L9");
    Grid global_grid2("L9");
    Grid g1(global_grid1, testdomain1);
    Grid g2(global_grid2, testdomain2);
    //Grid g1( global_grid1 );
    //Grid g2( global_grid2 );

    bool fourierTrc1 = true;
    bool fourierTrc2 = true;
    // fourierTrc1, fourierTrc2: need to be false if no global grid can be constructed
    //                           (like for grids created with LinearSpacing)

    //using LinearSpacing = grid::LinearSpacing;
    //StructuredGrid g2( LinearSpacing( {0., 10.}, 2 ), LinearSpacing( {-10., -90.}, 9 ) );
    // when using LinearSpacing: set fourierTrc2 to false

    int trc = 8;
    Trace t1(Here(), "translocal1 construction");
    trans::Trans transLocal1(global_grid1, g1.domain(), trc, option::type("local"));
    t1.stop();
    Trace t2(Here(), "translocal2 construction");
    //trans::Trans transLocal2( g2, trc, option::type( "local" ) );
    trans::Trans transLocal2(global_grid2, g2.domain(), trc, option::type("local"));
    t2.stop();

    double rav1 = 0., rav2 = 0.;  // compute average rms errors of transLocal1 and transLocal2

    functionspace::Spectral spectral(trc);

    int nb_scalar = 1, nb_vordiv = 1;
    int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
    std::vector<double> sp(2 * N * nb_scalar);
    std::vector<double> vor(2 * N * nb_vordiv);
    std::vector<double> div(2 * N * nb_vordiv);
    std::vector<double> rspecg(2 * N);
    std::vector<double> rgp1(nb_all * g1.size());
    std::vector<double> rgp2(nb_all * g2.size());
    std::vector<double> rgp1_analytic(g1.size());
    std::vector<double> rgp2_analytic(g2.size());

    StructuredGrid gs1(g1);
    StructuredGrid gs2(g2);
    double latPole   = 89.9999999;
    bool gridHasPole = gs1.y(0) > latPole || gs2.y(0) > latPole;
    double tolerance;
    int icase = 0;
    for (int ivar_in = 0; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
        for (int ivar_out = 0; ivar_out < 3; ivar_out++) {  // u, v, scalar
            int nb_fld = 1;
            if (ivar_out == 2 && !gridHasPole) {
                tolerance = 1.e-13;
                nb_fld    = nb_scalar;
            }
            else {
                tolerance = 2.e-6;
                nb_fld    = nb_vordiv;
            }
            for (int jfld = 0; jfld < nb_fld; jfld++) {  // multiple fields
                int k = 0;
                for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                    for (int n = m; n <= trc; n++) {             // total wavenumber
                        for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                            if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. &&
                                icase < 1000) {
                                auto start = std::chrono::system_clock::now();
                                for (int j = 0; j < 2 * N * nb_scalar; j++) {
                                    sp[j] = 0.;
                                }
                                for (int j = 0; j < 2 * N * nb_vordiv; j++) {
                                    vor[j] = 0.;
                                    div[j] = 0.;
                                }
                                if (ivar_in == 0) {
                                    vor[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 1) {
                                    div[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 2) {
                                    sp[k * nb_scalar + jfld] = 1.;
                                }

                                for (int j = 0; j < nb_all * g1.size(); j++) {
                                    rgp1[j] = 0.;
                                }
                                for (int j = 0; j < nb_all * g2.size(); j++) {
                                    rgp2[j] = 0.;
                                }
                                for (int j = 0; j < g1.size(); j++) {
                                    rgp1_analytic[j] = 0.;
                                }
                                for (int j = 0; j < g2.size(); j++) {
                                    rgp2_analytic[j] = 0.;
                                }

                                spectral_transform_grid_analytic(trc, fourierTrc1, n, m, imag, g1, rspecg.data(),
                                                                 rgp1_analytic.data(), ivar_in, ivar_out);

                                spectral_transform_grid_analytic(trc, fourierTrc2, n, m, imag, g2, rspecg.data(),
                                                                 rgp2_analytic.data(), ivar_in, ivar_out);

                                //Log::info() << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                //            << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;

                                ATLAS_TRACE_SCOPE("translocal1")
                                EXPECT_NO_THROW(transLocal1.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp1.data()));

                                int pos = (ivar_out * nb_vordiv + jfld);

                                double rms_gen1 =
                                    compute_rms(g1.size(), rgp1.data() + pos * g1.size(), rgp1_analytic.data());

                                //Log::info() << "rgp1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1[pos * g1.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic1:";
                                //for ( int j = 0; j < g1.size(); j++ ) {
                                //    Log::info() << rgp1_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl;

                                ATLAS_TRACE_SCOPE("translocal2")
                                EXPECT_NO_THROW(transLocal2.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp2.data()));

                                double rms_gen2 =
                                    compute_rms(g2.size(), rgp2.data() + pos * g2.size(), rgp2_analytic.data());

                                //Log::info() << std::endl << "rgp2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2[pos * g2.size() + j] << " ";
                                //};
                                //Log::info() << std::endl << "analytic2:";
                                //for ( int j = 0; j < g2.size(); j++ ) {
                                //    Log::info() << rgp2_analytic[j] << " ";
                                //};
                                //Log::info() << std::endl;
                                rav1 += rms_gen1;
                                rav2 += rms_gen2;
                                if (!(rms_gen1 < tolerance) || !(rms_gen2 < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
                                EXPECT(rms_gen1 < tolerance);
                                EXPECT(rms_gen2 < tolerance);
                                icase++;
                                auto end                                      = std::chrono::system_clock::now();  //
                                std::chrono::duration<double> elapsed_seconds = end - start;
                                std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                                std::string time_str = std::ctime(&end_time);
                                Log::debug() << "case " << icase << ", elapsed time: " << elapsed_seconds.count()
                                             << "s. Now: " << time_str.substr(0, time_str.length() - 1) << std::endl;
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }
    Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
    rav1 /= icase;
    Log::info() << "average RMS error of transLocal1: " << rav1 << std::endl;
    rav2 /= icase;
    Log::info() << "average RMS error of transLocal2: " << rav2 << std::endl;
}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE("test_trans_unstructured") {
    Log::info() << "test_trans_unstructured" << std::endl;
    // test transgeneral by comparing with analytic solution on an unstructured grid

    double tolerance = 1.e-13;

    //Domain testdomain = RectangularDomain( {20., 25.}, {40., 60.} );
    Domain testdomain = RectangularDomain({0., 90.}, {0., 90.});
    // Grid: (Adjust the following line if the test takes too long!)
    Grid grid_global("F32");
    Grid g(grid_global, testdomain);
    int trc = 31;
    StructuredGrid gs(g);
    std::vector<PointXY> pts(g.size());
    idx_t idx(0);
    for (idx_t j = 0; j < gs.ny(); ++j) {
        double lat = gs.y(j);
        for (idx_t i = 0; i < gs.nx(j); ++i) {
            double lon = gs.x(i, j);
            if (i == j && lat > 0) {
                //Log::info() << "idx=" << idx << " lon=" << lon << " lat=" << lat << std::endl;
                pts[idx++].assign(lon, lat);
            }
        }
    }
    Grid gu = UnstructuredGrid(new std::vector<PointXY>(&pts[0], &pts[idx]));
    Log::info() << "gu: size=" << gu.size() << std::endl;
    double rav1 = 0., rav2 = 0.;  // compute average rms errors of transLocal1 and transLocal2

    int nb_scalar = 1, nb_vordiv = 1;
    int N = (trc + 2) * (trc + 1) / 2, nb_all = nb_scalar + 2 * nb_vordiv;
    std::vector<double> sp(2 * N * nb_scalar);
    std::vector<double> vor(2 * N * nb_vordiv);
    std::vector<double> div(2 * N * nb_vordiv);
    std::vector<double> rspecg(2 * N);
    std::vector<double> gp(nb_all * g.size());
    std::vector<double> rgp1(nb_all * g.size());
    std::vector<double> rgp2(nb_all * g.size());
    std::vector<double> rgp_analytic1(g.size());
    std::vector<double> rgp_analytic2(gu.size());

    trans::Trans transLocal1(grid_global, testdomain, trc, option::type("local"));

    // ATLAS-173 : This should also work with precompute = true, ans should give same.
    trans::Trans transLocal2(gu, trc, option::type("local") | util::Config("precompute", false));

    int icase = 0;
    for (int ivar_in = 0; ivar_in < 3; ivar_in++) {         // vorticity, divergence, scalar
        for (int ivar_out = 0; ivar_out < 3; ivar_out++) {  // u, v, scalar
            int nb_fld = 1;
            if (ivar_out == 2) {
                tolerance = 1.e-13;
                nb_fld    = nb_scalar;
            }
            else {
                tolerance = 2.e-6;
                nb_fld    = nb_vordiv;
            }
            for (int jfld = 0; jfld < nb_fld; jfld++) {  // multiple fields
                int k = 0;
                for (int m = 0; m <= trc; m++) {                 // zonal wavenumber
                    for (int n = m; n <= trc; n++) {             // total wavenumber
                        for (int imag = 0; imag <= 1; imag++) {  // real and imaginary part

                            if (sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. &&
                                icase < 1000) {
                                auto start = std::chrono::system_clock::now();
                                for (int j = 0; j < 2 * N * nb_scalar; j++) {
                                    sp[j] = 0.;
                                }
                                for (int j = 0; j < 2 * N * nb_vordiv; j++) {
                                    vor[j] = 0.;
                                    div[j] = 0.;
                                }
                                if (ivar_in == 0) {
                                    vor[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 1) {
                                    div[k * nb_vordiv + jfld] = 1.;
                                }
                                if (ivar_in == 2) {
                                    sp[k * nb_scalar + jfld] = 1.;
                                }

                                for (int j = 0; j < nb_all * g.size(); j++) {
                                    gp[j]   = 0.;
                                    rgp1[j] = 0.;
                                    rgp2[j] = 0.;
                                }
                                for (int j = 0; j < g.size(); j++) {
                                    rgp_analytic1[j] = 0.;
                                }

                                for (int j = 0; j < gu.size(); j++) {
                                    rgp_analytic2[j] = 0.;
                                }

                                spectral_transform_grid_analytic(trc, false, n, m, imag, g, rspecg.data(),
                                                                 rgp_analytic1.data(), ivar_in, ivar_out);

                                //Log::info() << icase << " m=" << m << " n=" << n << " imag=" << imag << " structured: ";
                                ATLAS_TRACE_SCOPE("structured")
                                EXPECT_NO_THROW(transLocal1.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp1.data()));

                                spectral_transform_grid_analytic(trc, false, n, m, imag, gu, rspecg.data(),
                                                                 rgp_analytic2.data(), ivar_in, ivar_out);

                                //Log::info() << icase << " m=" << m << " n=" << n << " imag=" << imag << " unstructured: ";
                                ATLAS_TRACE_SCOPE("unstructured")
                                EXPECT_NO_THROW(transLocal2.invtrans(nb_scalar, sp.data(), nb_vordiv, vor.data(),
                                                                     div.data(), rgp2.data()));

                                int pos = (ivar_out * nb_vordiv + jfld);

                                double rms_gen1 =
                                    compute_rms(g.size(), rgp1.data() + pos * g.size(), rgp_analytic1.data());

                                double rms_gen2 =
                                    compute_rms(gu.size(), rgp2.data() + pos * gu.size(), rgp_analytic2.data());

                                rav1 += rms_gen1;
                                rav2 += rms_gen2;
                                if (!(rms_gen1 < tolerance) || !(rms_gen2 < tolerance)) {
                                    Log::info()
                                        << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out
                                        << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                    ATLAS_DEBUG_VAR(rms_gen1);
                                    ATLAS_DEBUG_VAR(rms_gen2);
                                    ATLAS_DEBUG_VAR(tolerance);
                                }
                                EXPECT(rms_gen1 < tolerance);
                                EXPECT(rms_gen2 < tolerance);
                                icase++;
                                auto end                                      = std::chrono::system_clock::now();  //
                                std::chrono::duration<double> elapsed_seconds = end - start;
                                std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                                std::string time_str = std::ctime(&end_time);
                                Log::debug() << "case " << icase << ", elapsed time: " << elapsed_seconds.count()
                                             << "s. Now: " << time_str.substr(0, time_str.length() - 1) << std::endl;
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }
    Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
    rav1 /= icase;
    Log::info() << "average RMS error of transLocal1: " << rav1 << std::endl;
    rav2 /= icase;
    Log::info() << "average RMS error of transLocal2: " << rav2 << std::endl;
}
#endif

//-----------------------------------------------------------------------------
#if ATLAS_HAVE_TRANS
CASE("test_trans_levels") {
    Log::info() << "test_trans_levels" << std::endl;

    // test trans_levels test puts
    //   the real component of spherical harmonic (m=1, n=2) on level 0
    //   and the imaginary component of spherical harmomic (m=2, n=3) on level 1
    //   of a regular Gaussian field with two levels.

    // it then runs directtrans and checks that the spectral coefficients are
    //   correct.
    //   (Note that this test is assuming we are running on 1PE.)

    // it finally runs the inverse transform and checks that we end up where we
    // started.

    std::string grid_uid("F2");  // Regular Gaussian F ( 8 N^2)
    StructuredGrid g(grid_uid);

    auto N             = GaussianGrid(g).N();  // -> cast to Gaussian grid and get the Gaussian number
    std::size_t levels = 2;
    std::vector<int> n{2, 3};     // total wavenumber on levels O and 1
    std::vector<int> m{1, 2};     // meridional wavenumber on levels O and 1
    std::vector<int> imag{0, 1};  // imaginary component or not
    const std::vector<idx_t> level_sequence{0, 1};

    functionspace::Spectral specFS(2 * N - 1, atlas::option::levels(levels));
    functionspace::StructuredColumns gridFS(g, atlas::option::levels(levels));

    atlas::trans::Trans transIFS(gridFS, specFS);

    Log::info() << "transIFS backend" << transIFS.backend() << std::endl;

    std::vector<atlas::PointLonLat> pointsLonLat;
    pointsLonLat.reserve(static_cast<std::size_t>(g.size()));
    for (auto ll : g.lonlat()) {
        pointsLonLat.push_back(ll);
    }
    atlas::Field gpf  = gridFS.createField<double>(atlas::option::name("gpf"));
    atlas::Field gpf2 = gridFS.createField<double>(atlas::option::name("gpf2"));
    atlas::Field spf  = specFS.createField<double>();

    auto view = atlas::array::make_view<double, 2>(gpf);
    for (atlas::idx_t j = gridFS.j_begin(); j < gridFS.j_end(); ++j) {
        for (atlas::idx_t i = gridFS.i_begin(j); i < gridFS.i_end(j); ++i) {
            atlas::idx_t jn = gridFS.index(i, j);
            for (atlas::idx_t jl : level_sequence) {
                view(jn, jl) =
                    sphericalharmonics_analytic_point(n[jl], m[jl], imag[jl], ((pointsLonLat[jn].lon() * M_PI) / 180.),
                                                      ((pointsLonLat[jn].lat() * M_PI) / 180.), 2, 2);
            }
        }
    }

    // transform fields to spectral and view
    transIFS.dirtrans(gpf, spf);

    auto spView = atlas::array::make_view<double, 2>(spf);

    Log::info() << "trans spectral coefficients = " << transIFS.spectralCoefficients() << std::endl;
    Log::info() << "checking the spectral coefficients are correct " << std::endl;

    int imagSize = 2;
    std::vector<int> index;
    for (int jl : level_sequence) {
        int offset(0);
        for (int k = 0; k < m[jl]; ++k) {
            offset += (2 * N - k) * imagSize;
        }
        index.push_back(offset + imag[jl] + (n[jl] - m[jl]) * imagSize);
    }

    for (int i = 0; i < transIFS.spectralCoefficients(); ++i) {
        for (int jl : level_sequence) {
            double value = (i == index[jl] ? 1.0 : 0.0);
            EXPECT(std::abs(spView(i, jl) - value) < 1e-14);
        }
    }

    transIFS.invtrans(spf, gpf2);

    Log::info() << "comparing going through direct and inverse transforms " << std::endl;

    auto view2 = atlas::array::make_view<double, 2>(gpf2);
    for (atlas::idx_t j = gridFS.j_begin(); j < gridFS.j_end(); ++j) {
        for (atlas::idx_t i = gridFS.i_begin(j); i < gridFS.i_end(j); ++i) {
            atlas::idx_t jn = gridFS.index(i, j);
            for (idx_t jl : level_sequence) {
                EXPECT(std::abs(view(jn, jl) - view2(jn, jl)) < 1e-14);
            }
        }
    }
}
#endif


#if ATLAS_HAVE_TRANS && (ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ))
CASE("test_2level_adjoint_test_with_powerspectrum_convolution") {
    std::string grid_uid("F64");  // Regular Gaussian F ( 8 N^2)
    atlas::StructuredGrid g2(grid_uid);

    auto N = atlas::GaussianGrid(g2).N();  // -> cast to Gaussian grid and get the Gaussian number

    auto levels = 2;
    std::vector<int> n{2, 3};            // total wavenumber on levels O and 1
    std::vector<int> m{1, 2};            // meridional wavenumber on levels O and 1
    std::vector<int> imag{false, true};  // imaginary component or not

    atlas::functionspace::Spectral specFS(2 * N - 1, atlas::option::levels(levels));

    atlas::functionspace::StructuredColumns gridFS(
        g2, atlas::grid::Partitioner(new atlas::grid::detail::partitioner::TransPartitioner()),
        atlas::option::levels(levels));

    auto lonlatview = atlas::array::make_view<double, 2>(gridFS.lonlat());

    atlas::trans::Trans transIFS(gridFS, specFS);

    //  Log::info() << "transIFS llatlon = "  << dynamic_cast<atlas::trans::TransIFS *>(transIFS.get())->trans()->llatlon << std::endl;

    Log::info() << "transIFS backend" << transIFS.backend() << std::endl;

    std::vector<float> powerSpectrum(2 * N, 0.0);

    for (std::size_t w = 0; w < powerSpectrum.size(); ++w) {
        powerSpectrum[w] = 1.0 / static_cast<float>(w + 1);
    }
    float tot = std::accumulate(powerSpectrum.cbegin(), powerSpectrum.cend(), 0.0);
    for (std::size_t w = 0; w < powerSpectrum.size(); ++w) {
        powerSpectrum[w] = powerSpectrum[w] / tot;
    }
    Log::info() << "create a fictitous power spectrum" << atlas::mpi::rank() << " " << powerSpectrum[0] << " "
                << powerSpectrum[1] << " " << tot << std::endl;

    atlas::Field gpf  = gridFS.createField<double>(atlas::option::name("gpf"));
    atlas::Field gpf2 = gridFS.createField<double>(atlas::option::name("gpf2"));
    atlas::Field spf  = specFS.createField<double>(atlas::option::name("spf"));
    atlas::Field spfg = specFS.createField<double>(atlas::option::name("spfg") | atlas::option::global());

    auto gpfView = atlas::array::make_view<double, 2>(gpf);
    for (atlas::idx_t j = gridFS.j_begin(); j < gridFS.j_end(); ++j) {
        for (atlas::idx_t i = gridFS.i_begin(j); i < gridFS.i_end(j); ++i) {
            atlas::idx_t jn = gridFS.index(i, j);
            for (atlas::idx_t jl : {0, 1}) {
                gpfView(jn, jl) = 0.0;
                if ((j == gridFS.j_end() - 1) && (i == gridFS.i_end(gridFS.j_end() - 1) - 1) && (jl == 0) &&
                    (atlas::mpi::rank() == 0)) {
                    gpfView(jn, jl) = 1.0;
                }
                if ((j == gridFS.j_end() - 1) && (i == gridFS.i_end(gridFS.j_end() - 1) - 1) && (jl == 1) &&
                    (atlas::mpi::rank() == 1)) {
                    gpfView(jn, jl) = 1.0;
                }
            }
        }
    }

    // transform fields to spectral and view
    transIFS.invtrans_adj(gpf, spf);

    auto spfView                   = atlas::array::make_view<double, 2>(spf);
    const auto zonal_wavenumbers   = specFS.zonal_wavenumbers();
    const int nb_zonal_wavenumbers = zonal_wavenumbers.size();
    double adj_value(0.0);
    int i{0};
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
        const std::size_t m1 = zonal_wavenumbers(jm);
        for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(2 * N - 1); ++n1) {
            for (int imag1 : {0, 1}) {
                for (int jl : {0, 1}) {
                    // scale by the square root of the power spectrum
                    // note here that we need the modal power spectrum
                    // i.e. divide powerSpectra by (2n +1)
                    spfView(i, jl) = spfView(i, jl) * std::sqrt(std::sqrt(static_cast<double>(powerSpectrum[n1]))) /
                                     std::sqrt(static_cast<double>(2 * n1 + 1));

                    // adjoint at the heart
                    double temp = spfView(i, jl) * spfView(i, jl);
                    adj_value += (m1 > 0 ? 2 * temp : temp);  // to take account of -ve m.

                    // scale by the square root of the power spectrum (again!)
                    spfView(i, jl) = spfView(i, jl) * std::sqrt(std::sqrt(static_cast<double>(powerSpectrum[n1]))) /
                                     std::sqrt(static_cast<double>(2 * n1 + 1));
                }
                ++i;
            }
        }
    }
    atlas::mpi::comm().allReduceInPlace(adj_value, eckit::mpi::sum());

    transIFS.invtrans(spf, gpf2);

    Log::info() << "adjoint test transforms " << std::endl;
    auto gpf2View = atlas::array::make_view<double, 2>(gpf2);
    double adj_value2(0.0);
    for (atlas::idx_t j = gridFS.j_begin(); j < gridFS.j_end(); ++j) {
        for (atlas::idx_t i = gridFS.i_begin(j); i < gridFS.i_end(j); ++i) {
            atlas::idx_t jn = gridFS.index(i, j);
            for (atlas::idx_t jl : {0, 1}) {
                adj_value2 += gpfView(jn, jl) * gpf2View(jn, jl);
            }
        }
    }
    atlas::mpi::comm().allReduceInPlace(adj_value2, eckit::mpi::sum());

    Log::info() << "adjoint test "
                << " " << adj_value << " " << adj_value2 << std::endl;
    EXPECT(std::abs(adj_value - adj_value2) / adj_value < 1e-12);
}
#endif


#if 0
CASE( "test_trans_fourier_truncation" ) {
    Log::info() << "test_trans_fourier_truncation" << std::endl;
    // test transgeneral by comparing its result with the trans library
    // this test is based on the test_nomesh case in test_trans.cc

    Grid g( "F640" );
    grid::StructuredGrid gs( g );
    int ndgl = gs.ny();
    //int trc = 2*ndgl; // extreme high truncation (below linear)
    int trc = ndgl - 1;  // linear
    //int trc = 5./6.*ndgl-1; // between linear and quadratic
    //int trc = 2./3.*ndgl-1; // quadratic
    //int trc = ndgl/2. -1; // cubic
    trans::Trans trans( g, trc );
    for ( int j = 0; j < gs.ny(); j += 80 ) {
        double lat = gs.y( j ) * util::Constants::degreesToRadians();
        int trcFT  = trans::fourier_truncation( trc, gs.nx( j ), gs.nxmax(), gs.ny(), lat, grid::RegularGrid( g ) );
        Log::info() << trcFT << "         " << gs.nx( j ) << std::endl;
    }
    // TODO: create some real criterion to test fourier_truncation. So far only comparison with trans library through print statements.
}
#endif
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run<atlas::test::AtlasTransEnvironment>(argc, argv);
}
