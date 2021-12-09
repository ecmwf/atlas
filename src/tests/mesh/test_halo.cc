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

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/library/config.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/Unique.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/TestMeshes.h"

using namespace atlas::output;
using namespace atlas::util;
using namespace atlas::meshgenerator;

namespace atlas {
namespace test {

double dual_volume(Mesh& mesh) {
    mesh::Nodes& nodes = mesh.nodes();
    mesh::IsGhostNode is_ghost_node(nodes);
    int nb_nodes                             = nodes.size();
    array::ArrayView<double, 1> dual_volumes = array::make_view<double, 1>(nodes.field("dual_volumes"));
    double area                              = 0;
    for (int node = 0; node < nb_nodes; ++node) {
        if (!is_ghost_node(node)) {
            area += dual_volumes(node);
        }
    }

    ATLAS_TRACE_MPI(ALLREDUCE) { mpi::comm().allReduceInPlace(area, eckit::mpi::sum()); }

    return area;
}

#if 0
CASE( "test_small" )
{
  int nlat = 5;
  int lon[5] = {10, 12, 14, 16, 16};

  Mesh = test::generate_mesh(nlat, lon);

  mesh::actions::build_parallel_fields(*m);
  mesh::actions::build_periodic_boundaries(*m);
  mesh::actions::build_halo(*m,2);


  if( mpi::comm().size() == 5 )
  {
    IndexView<int,1> ridx ( m->nodes().remote_index() );
    array::array::ArrayView<gidx_t,1> gidx ( m->nodes().global_index() );

    switch( mpi::comm().rank() ) // with 5 tasks
    {
    case 0:
      EXPECT( ridx(9) ==  9  );
      EXPECT( gidx(9) ==  10 );
      EXPECT( ridx(30) == 9 );
      EXPECT( gidx(30) == 875430066 ); // hashed unique idx
      break;
    }
  }
  else
  {
    if( mpi::comm().rank() == 0 )
      std::cout << "skipping tests with 5 mpi tasks!" << std::endl;
  }

  mesh::actions::build_edges(*m);
  mesh::actions::build_median_dual_mesh(*m);

  EXPECT( eckit::types::is_approximately_equal( test::dual_volume(*m), 2.*M_PI*M_PI, 1e-6 ));

  std::stringstream filename; filename << "small_halo_p" << mpi::comm().rank() << ".msh";
  Gmsh(filename.str()).write(*m);
}
#endif

#if 1
CASE("test_custom") {
    // Mesh m = test::generate_mesh( T63() );

    Mesh m = test::generate_mesh({10, 12, 14, 16, 16, 16, 16, 14, 12, 10});

    mesh::actions::build_nodes_parallel_fields(m.nodes());
    mesh::actions::build_periodic_boundaries(m);
    mesh::actions::build_halo(m, 1);

    std::stringstream filename;
    filename << "custom.msh";
    Gmsh(filename.str(), util::Config("ghost", true)).write(m);

    //  EXPECT( eckit::types::is_approximately_equal( test::dual_volume(m),
    //  2.*M_PI*M_PI, 1e-6 ));

    auto lonlat = array::make_view<double, 2>(m.nodes().lonlat());

#if ATLAS_BITS_GLOBAL == 64

    std::vector<uidx_t> check;
    switch (mpi::comm().rank()) {
        case 0:
            check = {607990293346953216, 607990293382953216, 607990293418953216, 607990293454953216, 607990293490953216,
                     607990293526953216, 607990293562953216, 607990293598953216, 607990293634953216, 607990293670953216,
                     644481443595331584, 644481443625331584, 644481443655331584, 644481443685331584, 644481443715331584,
                     644481443745331584, 644481443775331584, 644481443805331584, 644481443835331584, 644481443865331584,
                     644481443895331584, 644481443925331584, 681187136050079744, 681187136075794030, 681187136101508315,
                     681187136127222601, 681187136152936887, 681187136178651173, 607990293706953216, 644481443955331584,
                     681187136204365458, 681187136230079744, 681187136255794030, 681187136281508315, 681187136307222601,
                     681187136332936887, 681187136358651173, 681187136384365458, 717939789677242368, 717939789699742368,
                     717939789722242368, 717939789744742368, 717939789767242368, 717939789789742368, 717939789812242368,
                     754708008265885696, 754708008288385696, 754708008310885696, 754708008333385696, 754708008355885696,
                     754708008378385696, 754708008400885696, 717939789834742368, 717939789857242368, 717939789879742368,
                     717939789902242368, 754708008423385696, 681187136410079744, 717939789924742368, 717939789947242368,
                     717939789969742368, 717939789992242368, 717939790014742368, 607990293310953216, 644481443565331584,
                     681187136024365458, 717939789654742368, 754708008243385696, 607990293742953216, 644481443985331584,
                     681187136435794030};
            break;
        case 1:
            check = {717939789677242368, 717939789699742368, 717939789722242368, 717939789744742368, 717939789767242368,
                     717939789789742368, 754708008265885696, 754708008288385696, 754708008310885696, 754708008333385696,
                     754708008355885696, 754708008378385696, 791480221174114304, 791480221196614304, 791480221219114304,
                     791480221241614304, 791480221264114304, 828248439762757632, 828248439785257632, 828248439807757632,
                     828248439830257632, 828248439852757632, 865001093389920256, 865001093415634542, 865001093441348827,
                     865001093467063113, 865001093492777399, 717939789812242368, 754708008400885696, 791480221286614304,
                     828248439875257632, 865001093518491685, 901706785844668416, 901706785874668416, 901706785904668416,
                     901706785934668416, 901706785964668416, 681187136050079744, 681187136075794030, 681187136101508315,
                     681187136127222601, 681187136152936887, 681187136178651173, 681187136204365458, 717939789834742368,
                     754708008423385696, 791480221309114304, 791480221331614304, 828248439897757632, 865001093544205970,
                     901706785994668416, 938197936093046784, 938197936129046784, 938197936165046784, 938197936201046784,
                     938197936237046784, 717939789654742368, 754708008243385696, 791480221151614304, 828248439740257632,
                     865001093364205970, 901706785814668416};
            break;
        case 2:
            check = {681187136204365458, 681187136230079744, 681187136255794030, 717939789812242368, 717939789834742368,
                     717939789857242368, 717939789879742368, 717939789902242368, 754708008400885696, 754708008423385696,
                     754708008445885696, 754708008468385696, 754708008490885696, 791480221286614304, 791480221309114304,
                     791480221331614304, 791480221354114304, 791480221376614304, 828248439875257632, 828248439897757632,
                     828248439920257632, 828248439942757632, 828248439965257632, 865001093518491685, 865001093544205970,
                     865001093569920256, 865001093595634542, 681187136281508315, 717939789924742368, 754708008378385696,
                     754708008513385696, 791480221399114304, 828248439987757632, 865001093492777399, 865001093621348827,
                     901706785964668416, 901706785994668416, 901706786024668416, 901706786054668416, 644481443745331584,
                     644481443775331584, 644481443805331584, 644481443835331584, 681187136178651173, 681187136307222601,
                     717939789789742368, 717939789767242368, 754708008355885696, 791480221264114304, 828248439852757632,
                     865001093467063113, 901706785934668416, 717939789947242368, 754708008535885696, 791480221421614304,
                     791480221444114304, 828248440010257632, 865001093647063113, 901706786084668416, 938197936201046784,
                     938197936237046784, 938197936273046784, 938197936309046784};
            break;
        case 3:
            check = {681187136281508315, 681187136307222601, 681187136332936887, 681187136358651173, 681187136384365458,
                     717939789924742368, 717939789947242368, 717939789969742368, 717939789992242368, 717939790014742368,
                     754708008513385696, 754708008535885696, 754708008558385696, 754708008580885696, 754708008603385696,
                     791480221399114304, 791480221421614304, 791480221444114304, 791480221466614304, 791480221489114304,
                     791480221511614304, 828248439987757632, 828248440010257632, 828248440032757632, 828248440055257632,
                     828248440077757632, 828248440100257632, 644481443955331584, 681187136410079744, 717939789902242368,
                     717939790037242368, 754708008490885696, 754708008625885696, 791480221534114304, 828248440122757632,
                     865001093621348827, 865001093647063113, 865001093672777399, 865001093698491685, 865001093724205970,
                     865001093749920256, 607990293706953216, 644481443805331584, 644481443835331584, 644481443865331584,
                     644481443895331584, 644481443925331584, 681187136255794030, 717939789879742368, 754708008468385696,
                     791480221376614304, 828248439965257632, 865001093595634542, 901706786054668416, 901706786084668416,
                     901706786114668416, 901706786144668416, 901706786174668416, 901706786204668416, 644481443985331584,
                     681187136435794030, 717939790059742368, 754708008648385696, 791480221556614304, 828248440145257632,
                     865001093775634542};
            break;
        case 4:
            check = {865001093621348827, 865001093647063113, 865001093672777399, 865001093698491685, 865001093724205970,
                     901706785844668416, 901706785874668416, 901706785904668416, 901706785934668416, 901706785964668416,
                     901706785994668416, 901706786024668416, 901706786054668416, 901706786084668416, 901706786114668416,
                     901706786144668416, 901706786174668416, 938197936093046784, 938197936129046784, 938197936165046784,
                     938197936201046784, 938197936237046784, 938197936273046784, 938197936309046784, 938197936345046784,
                     938197936381046784, 938197936417046784, 865001093749920256, 901706786204668416, 938197936453046784,
                     865001093389920256, 865001093415634542, 865001093441348827, 865001093467063113, 865001093492777399,
                     865001093518491685, 828248439987757632, 865001093544205970, 865001093569920256, 865001093595634542,
                     828248440010257632, 828248440032757632, 828248440055257632, 828248440077757632, 828248440100257632,
                     828248440122757632, 865001093364205970, 901706785814668416, 938197936057046784, 828248440145257632,
                     865001093775634542, 901706786234668416, 938197936489046784};
            break;
        default:
            check.clear();
    }
    std::vector<uidx_t> uid(m.nodes().size());
    for (idx_t j = 0; j < m.nodes().size(); ++j) {
        uid[j] = util::unique_lonlat(lonlat(j, 0), lonlat(j, 1));
    }
    if (check.size() && mpi::comm().size() == 5) {
        EXPECT_EQ(uid.size(), check.size());

        if (uid != check) {
            for (idx_t i = 0; i < uid.size(); ++i) {
                if (uid[i] != check[i]) {
                    Log::warning() << "uid[" << i << "] != check[" << i << "] : " << uid[i] << " expected to be "
                                   << check[i] << "   point = " << std::setprecision(7) << std::fixed
                                   << PointLonLat{lonlat(i, LON), lonlat(i, LAT)}
                                   << "   microdeg(lon) = " << microdeg(lonlat(i, LON)) << std::endl;
                }
            }
            Log::info() << "uid = { ";
            for (idx_t i = 0; i < uid.size(); ++i) {
                Log::info() << uid[i] << std::string(i < uid.size() - 1 ? ", " : " ");
            }
            Log::info() << "};" << std::endl;
        }
        EXPECT(uid == check);
    }
#endif
    //  FunctionSpace& edges = m.function_space("edges");
    //  array::array::ArrayView<double,1> dual_volumes  ( nodes.field(
    //  "dual_volumes" ) );
    //  array::array::ArrayView<double,2> dual_normals  ( edges.field(
    //  "dual_normals" ) );

    //  std::string checksum;
    //  checksum = nodes.checksum()->execute(dual_volumes);
    //  DEBUG("dual_volumes checksum "<<checksum,0);

    //  checksum = edges.checksum()->execute(dual_normals);
    //  DEBUG("dual_normals checksum "<<checksum,0);
}
#endif
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
