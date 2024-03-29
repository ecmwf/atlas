/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


//-----------------------------------------------------------------

#include <cmath>

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/Reorder.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"


#include "eckit/config/Resource.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/linalg/Triplet.h"

#include "tests/AtlasTestEnvironment.h"


#include "atlas/output/detail/GmshIO.h"

//-----------------------------------------------------------------

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

void outputConnectedCoordinates(const eckit::PathName& filepath, const array::ArrayView<double, 2>& xy,
                                const std::string& title) {
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    int mpi_rank                 = int(comm.rank());
    int mpi_size                 = int(comm.size());

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for (idx_t i = 0; i < xy.shape(0); ++i) {
        xmin = std::min(xmin, xy(i, XX));
        xmax = std::max(xmax, xy(i, XX));
    }
    comm.allReduceInPlace(xmin, eckit::mpi::min());
    comm.allReduceInPlace(xmax, eckit::mpi::max());

    idx_t count     = xy.shape(0);
    idx_t count_all = count;
    comm.allReduceInPlace(count_all, eckit::mpi::sum());

    for (int r = 0; r < mpi_size; ++r) {
        if (mpi_rank == r) {
            std::ofstream f(filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app);
            // clang-format off
            if ( mpi_rank == 0 ) {
                f << "\n"
                     "\n" "import matplotlib.pyplot as plt"
                     "\n" "from matplotlib.path import Path"
                     "\n" "import matplotlib.patches as patches"
                     "\n" ""
                     "\n" "from itertools import cycle"
                     "\n" "import matplotlib.cm as cm"
                     "\n" "import numpy as np"
                     "\n" "cycol = cycle([cm.Paired(i) for i in np.linspace(0,1,12,endpoint=True)]).next"
                     "\n" ""
                     "\n" "fig = plt.figure()"
                     "\n" "ax = fig.add_subplot(111,aspect='equal')"
                     "\n" "";
            }

            if ( mpi_rank == r ) {  // replace "r" with rank you wish to plot only
                f << "\n" "verts_"
                  << r << " = [";
                for ( idx_t n = 0; n < xy.shape( 0 ); ++n ) {
                    idx_t i = n;  //reorder.hilbert_reordering_[n].second;
                    f << "\n  (" << xy( i, XX ) << ", " << xy( i, YY ) << "), ";
                }
                f << "\n]"
                     "\n" ""
                     "\n" "codes_" << r << " = [Path.MOVETO]"
                     "\n" "codes_" << r << ".extend([Path.LINETO] * " << ( xy.shape( 0 ) - 2 ) << ")"
                     "\n" "codes_" << r << ".extend([Path.LINETO])"
                     "\n" ""
                     "\n" "count_" << r << " = " << count
                  << "\n" "count_all_" << r << " = " << count_all
                  << "\n" "c = cycol()"
                     "\n" "xs_" << r << ", ys_" << r << " = zip(*verts_" << r << ")"
                     "\n" "ax.plot(xs_" << r << ",ys_" << r << ", '-', lw=1, color=c )"
                     //"\n" "for i in range( len(verts_0) ):"
                     //"\n" "  plt.text( xs_0[i], ys_0[i], str(i) )"
                     "\n" "";
            }
            if ( mpi_rank == mpi_size - 1 ) {
                f << "\n" "ax.set_xlim( " << xmin << "-5, " << xmax << "+5)"
                     "\n" "ax.set_ylim(-90-5,  90+5)"
                     "\n" "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
                     "\n" "ax.set_yticks([-90,-45,0,45,90])"
                     "\n" "plt.grid()"
                     "\n" "plt.title('" << title << "')"
                     "\n" "plt.show()";
            }
            // clang-format on
        }
        comm.barrier();
    }
}

void outputTriplets(const eckit::PathName& filepath, const std::vector<eckit::linalg::Triplet>& triplets, idx_t rows,
                    idx_t cols, const std::string& title, const std::string& xlabel, const std::string& ylabel) {
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    int mpi_rank                 = int(comm.rank());
    int mpi_size                 = int(comm.size());

    for (int r = 0; r < mpi_size; ++r) {
        if (mpi_rank == r) {
            std::ofstream f(filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app);

            if (mpi_rank == 0) {
                // clang-format off
                f << "\n"
                     "\n"   "import matplotlib.pyplot as plt"
                     "\n"   "from matplotlib.path import Path"
                     "\n"   "import matplotlib.patches as patches"
                     "\n"   "import scipy.sparse as sparse"
                     "\n"   ""
                     "\n"   "from itertools import cycle"
                     "\n"   "import matplotlib.cm as cm"
                     "\n"   "import numpy as np"
                     "\n"   "cycol = cycle([cm.Paired(i) for i in np.linspace(0,1,12,endpoint=True)]).next"
                     "\n"   ""
                     "\n"   "fig = plt.figure()"
                     //"\n"   "ax = fig.add_subplot(111,aspect='equal')"
                     "\n"   "";
            }

            if ( mpi_rank == 0 ) {  // replace "r" with rank you wish to plot only
                f << "\n"   "row_"<< r << " = np.array([";
                for ( const auto& triplet : triplets ) {
                    f << "\n  " << triplet.row() << ", ";
                }
                f << "\n"   "])"
                     "\n"   ""
                     "\n"   "col_"<< r << " = np.array([";
                for ( const auto& triplet : triplets ) {
                    f << "\n  " << triplet.col() << ", ";
                }
                f << "\n"   "])"
                     "\n"   "data_"<< r << " = np.array([";
                for ( const auto& triplet : triplets ) {
                    f << "\n  " << triplet.value() << ", ";
                }
                f << "\n"   "])"
                     "\n"   ""
                     "\n"   "matrix_" << r << " = sparse.csr_matrix((data_"<<r<<",(row_"<<r<<",col_"<<r<<")),shape=("<<rows<<","<<cols<<"))"
                     "\n"   ""
                     "\n"   "plt.spy( matrix_"<<r<<",markersize=3 )"
                     "\n"   "";
            }
            if ( mpi_rank == mpi_size - 1 ) {
                f << "\n"   "plt.grid()"
                     "\n"   "plt.axes().set_aspect('auto', 'datalim')"
                     "\n"   "plt.title('" << title << "')"
                     "\n"   "plt.xlabel('" << xlabel << "')"
                     "\n"   "plt.ylabel('" << ylabel << "')"
                     "\n"   "plt.show()";
            }
            // clang-format on
        }
        comm.barrier();
    }
}

Field create_cell_centres(Mesh& mesh) {
    auto cell_centres =
        Field("cell_centres", array::make_datatype<double>(), array::make_shape(mesh.cells().size(), 2));
    auto nodes_xy = array::make_view<double, 2>(mesh.nodes().xy());
    for (idx_t t = 0; t < mesh.cells().nb_types(); ++t) {
        auto& cells = mesh.cells().elements(t);
        auto xy     = cells.view<double, 2>(cell_centres);

        // Compute cell-centres
        {
            const auto& node_connectivity = cells.node_connectivity();
            const idx_t nb_nodes          = cells.nb_nodes();
            const double nb_nodes_double  = nb_nodes;
            for (idx_t e = 0; e < cells.size(); ++e) {
                double x{0};
                double y{0};
                for (idx_t c = 0; c < nb_nodes; ++c) {
                    idx_t n = node_connectivity(e, c);
                    x += nodes_xy(n, XX);
                    y += nodes_xy(n, YY);
                }
                xy(e, XX) = x / nb_nodes_double;
                xy(e, YY) = y / nb_nodes_double;
            }
        }
    }
    return cell_centres;
}

Field create_edge_centres(Mesh& mesh) {
    auto edge_centres =
        Field("edge_centres", array::make_datatype<double>(), array::make_shape(mesh.edges().size(), 2));
    auto nodes_xy = array::make_view<double, 2>(mesh.nodes().xy());
    auto& edges   = mesh.edges();
    auto xy       = array::make_view<double, 2>(edge_centres);

    // Compute edge-centres
    {
        const auto& node_connectivity = edges.node_connectivity();
        for (idx_t e = 0; e < edges.size(); ++e) {
            double x{0};
            double y{0};
            for (idx_t c = 0; c < 2; ++c) {
                idx_t n = node_connectivity(e, c);
                x += nodes_xy(n, XX);
                y += nodes_xy(n, YY);
            }
            xy(e, XX) = 0.5 * x;
            xy(e, YY) = 0.5 * y;
        }
    }
    return edge_centres;
}

std::string grid_name() {
    return eckit::Resource<std::string>("--grid", "O16");
}
Mesh get_mesh() {
    if (grid_name() != "unstructured") {
        auto generate_mesh = StructuredMeshGenerator(util::Config("patch_pole", false)("triangulate", true));
        return generate_mesh(Grid{grid_name()});
    }
    else {
        output::detail::GmshIO gmsh_reader;
        std::string file = eckit::Resource<std::string>("--mesh", "");
        return gmsh_reader.read(file);
    }
}

void test_reordering(const util::Config& reorder_config, const std::vector<idx_t>& expected = {}) {
    std::string type = reorder_config.getString("type");
    // clang-format off
    std::string title = type == "none" ? "Default order" :
                        type == "hilbert" ? "Hilbert space filling curve" :
                        type == "reverse_cuthill_mckee" ? "Reverse Cuthill Mckee order":
                        type;
    type = grid_name() + "_" + type;
    // clang-format on
    auto mesh = get_mesh();

    auto reorder = mesh::actions::Reorder{reorder_config};

    if (expected.size()) {
        const auto node_order = reorder.get()->computeNodesOrder(mesh);
        EXPECT(node_order == expected);
    }

    if (false) {  // output node_order to feed in to expected
        const auto _node_order = reorder.get()->computeNodesOrder(mesh);
        Log::info() << "{";
        for (size_t i = 0; i < _node_order.size(); ++i) {
            Log::info() << _node_order[i] << ",";
        }
        Log::info() << "}";
    }

    reorder(mesh);

    bool sort_edges = false;
    //    bool sort_edges = false;
    mesh::actions::build_edges(mesh, util::Config("sort_edges", sort_edges));

    auto xy = array::make_view<double, 2>(mesh.nodes().xy());
    outputConnectedCoordinates(type + "_nodes.py", xy, title);
    output::Gmsh gmsh{type + ".msh", util::Config("coordinates", "xy")};
    gmsh.write(mesh);

    Field cell_centres   = create_cell_centres(mesh);
    auto cell_centres_xy = array::make_view<double, 2>(cell_centres);
    outputConnectedCoordinates(type + "_elements.py", cell_centres_xy, title);

    Field edge_centres   = create_edge_centres(mesh);
    auto edge_centres_xy = array::make_view<double, 2>(edge_centres);
    outputConnectedCoordinates(type + "_edges.py", edge_centres_xy, title);

    std::vector<eckit::linalg::Triplet> n2n_triplets;
    std::vector<eckit::linalg::Triplet> e2n_triplets;
    const auto& node_connectivity = mesh.edges().node_connectivity();
    const idx_t nb_edges          = mesh.edges().size();
    n2n_triplets.reserve(2 * nb_edges);
    for (idx_t e = 0; e < nb_edges; ++e) {
        n2n_triplets.emplace_back(node_connectivity(e, 0), node_connectivity(e, 1), 1.);
        n2n_triplets.emplace_back(node_connectivity(e, 1), node_connectivity(e, 0), 1.);
        e2n_triplets.emplace_back(e, node_connectivity(e, 0), 1.);
        e2n_triplets.emplace_back(e, node_connectivity(e, 1), 1.);
    }
    std::sort(n2n_triplets.begin(), n2n_triplets.end());
    outputTriplets(type + "_n2n_triplets.py", n2n_triplets, mesh.nodes().size(), mesh.nodes().size(), title, "nodes",
                   "nodes");
    outputTriplets(type + "_e2n_triplets.py", e2n_triplets, mesh.edges().size(), mesh.nodes().size(), title,
                   "connected nodes", "edges");
}

CASE("test_hilbert_reordering") {
    auto reorder_config = option::type("hilbert") | util::Config("recursion", 30);

    std::vector<idx_t> expected;
    if (grid_name() == "O16" && mpi::comm().size() == 1) {
        // clang-format off
        expected = {0,20,1,21,45,74,73,44,72,104,105,141,140,180,224,225,181,182,226,227,142,106,107,143,183,228,184,229,230,185,145,144,108,76,47,75,46,22,2,23,3,77,48,78,49,24,4,25,5,26,51,80,79,50,111,148,112,113,150,149,190,191,236,235,234,189,233,232,187,188,147,110,146,109,186,231,279,280,332,331,388,448,449,450,389,390,451,452,391,334,333,281,282,335,336,283,284,285,338,337,394,395,455,456,454,453,393,392,517,586,587,518,519,520,521,589,590,588,661,662,663,739,740,738,737,660,659,735,736,734,733,656,657,658,585,516,515,584,583,514,513,581,582,654,655,732,731,730,729,653,652,728,727,725,726,649,650,651,578,577,509,510,511,579,580,512,447,446,386,387,278,277,330,329,276,328,384,385,445,444,443,442,382,383,327,275,274,326,325,273,272,324,380,381,441,440,504,505,572,573,574,506,507,508,576,575,647,648,724,723,722,646,645,644,720,721,800,801,802,882,881,880,957,956,1028,1029,1030,958,959,1031,1032,960,884,883,803,804,805,806,885,886,887,807,808,810,809,889,888,963,965,964,1036,1035,1034,962,961,1033,1101,1100,1164,1165,1102,1103,1104,1167,1166,1226,1227,1282,1281,1225,1224,1280,1279,1278,1222,1223,1163,1099,1098,1162,1161,1097,1096,1160,1220,1221,1277,1276,1328,1376,1377,1329,1330,1331,1378,1422,1462,1461,1421,1420,1460,1496,1497,1528,1556,1580,1557,1581,1529,1498,1499,1530,1558,1582,1559,1531,1500,1464,1425,1424,1463,1423,1380,1379,1332,1333,1381,1334,1382,1335,1336,1337,1384,1383,1427,1466,1426,1465,1501,1532,1583,1560,1533,1502,1503,1534,1584,1561,1585,1562,1535,1504,1469,1430,1429,1468,1467,1428,1385,1338,1339,1386,1387,1340,1289,1288,1233,1234,1175,1112,1111,1110,1174,1173,1109,1108,1172,1231,1232,1287,1286,1285,1230,1229,1284,1283,1228,1169,1168,1105,1106,1107,1170,1171,1040,969,968,1039,1038,1037,966,967,891,890,811,812,813,892,893,894,815,814,816,817,896,895,970,1041,1042,971,972,1043,1044,1045,974,973,899,898,897,818,819,820,821,822,901,900,975,1046,1047,976,977,1048,1049,978,903,902,823,824,825,826,904,905,906,827,828,830,829,908,907,981,983,982,1053,1052,1051,980,979,1050,1117,1116,1179,1180,1118,1119,1120,1182,1181,1240,1241,1295,1294,1239,1238,1293,1292,1291,1236,1237,1178,1115,1114,1177,1176,1113,1235,1290,1388,1341,1342,1343,1389,1432,1471,1470,1431,1505,1563,1586,1536,1506,1507,1537,1564,1587,1565,1538,1508,1473,1435,1434,1472,1433,1391,1390,1344,1345,1392,1346,1393,1347,1348,1349,1395,1394,1437,1475,1436,1474,1509,1539,1588,1566,1540,1510,1511,1541,1589,1567,1439,1477,1476,1438,1396,1350,1351,1397,1301,1247,1127,1126,1189,1188,1125,1124,1187,1245,1246,1300,1299,1298,1244,1243,1297,1296,1242,1184,1183,1121,1122,1123,1185,1186,1057,987,986,1056,1055,1054,984,985,910,909,831,832,833,911,912,913,835,834,836,837,915,914,988,1058,1059,989,990,1060,1061,991,917,916,838,839,759,681,680,758,757,755,756,678,679,605,604,534,535,536,606,607,537,471,470,409,297,351,350,296,295,349,407,408,469,468,467,466,405,406,348,294,347,346,293,292,345,403,404,465,464,530,599,600,601,531,532,533,603,602,675,676,677,754,753,752,674,673,750,751,749,748,671,672,598,529,528,597,596,527,526,595,668,669,670,747,746,745,744,667,666,743,742,741,664,665,592,591,522,523,524,593,594,525,459,460,398,397,458,457,396,339,286,287,340,341,288,289,343,342,399,461,462,400,401,463,402,344,290,291,242,241,196,155,117,154,116,153,194,195,240,239,193,238,237,192,151,114,152,115,82,52,81,27,6,28,53,83,54,84,29,7,8,30,56,86,85,55,118,156,197,243,198,244,245,199,157,119,120,158,246,200,201,247,159,121,87,57,31,9,10,32,11,33,59,90,89,58,88,122,123,161,160,202,248,249,203,204,250,251,162,124,125,163,205,252,206,253,254,207,165,164,126,92,61,91,60,34,12,35,13,93,62,94,63,36,14,37,15,38,65,96,95,64,129,168,130,131,170,169,212,213,260,259,258,211,257,256,209,210,167,128,166,127,208,255,305,306,360,359,418,480,481,482,419,420,483,484,421,362,361,307,308,363,364,309,310,311,366,365,424,425,487,488,486,485,423,422,551,622,623,552,553,554,555,625,626,624,699,700,701,779,780,778,777,698,697,775,776,774,773,694,695,696,621,550,549,620,619,548,547,617,618,692,693,772,771,770,769,691,690,768,767,765,766,687,688,689,614,613,543,544,545,615,616,546,479,478,416,417,304,303,358,357,302,356,414,415,477,476,475,474,412,413,355,301,300,354,353,299,298,352,410,411,473,472,538,539,608,609,610,540,541,542,612,611,685,686,764,763,762,684,683,682,760,761,840,841,842,920,919,918,993,992,1062,1063,1064,994,995,1065,1066,996,922,921,843,844,845,846,923,924,925,847,848,850,849,927,926,999,1001,1000,1070,1069,1068,998,997,1067,1133,1132,1194,1195,1134,1135,1136,1197,1196,1254,1255,1308,1307,1253,1252,1306,1305,1304,1250,1251,1193,1131,1130,1192,1191,1129,1128,1190,1248,1249,1303,1302,1352,1398,1399,1353,1354,1355,1400,1442,1480,1479,1441,1440,1478,1512,1513,1542,1568,1590,1569,1591,1543,1514,1515,1544,1570,1592,1571,1545,1516,1482,1445,1444,1481,1443,1402,1401,1356,1357,1403,1358,1404,1359,1360,1361,1406,1405,1447,1484,1446,1483,1517,1546,1593,1572,1547,1518,1519,1548,1594,1573,1595,1574,1549,1520,1487,1450,1449,1486,1485,1448,1407,1362,1363,1408,1409,1364,1315,1314,1261,1262,1205,1144,1143,1142,1204,1203,1141,1140,1202,1259,1260,1313,1312,1311,1258,1257,1310,1309,1256,1199,1198,1137,1138,1139,1200,1201,1074,1005,1004,1073,1072,1071,1002,1003,929,928,851,852,853,930,931,932,855,854,856,857,934,933,1006,1075,1076,1007,1008,1077,1078,1079,1010,1009,937,936,935,858,859,860,861,862,939,938,1011,1080,1081,1012,1013,1082,1083,1014,941,940,863,864,865,866,942,943,944,867,868,870,869,946,945,1017,1019,1018,1087,1086,1085,1016,1015,1084,1149,1148,1209,1210,1150,1151,1152,1212,1211,1268,1269,1321,1320,1267,1266,1319,1318,1317,1264,1265,1208,1147,1146,1207,1206,1145,1263,1316,1410,1365,1366,1367,1411,1452,1489,1488,1451,1521,1575,1596,1550,1522,1523,1551,1576,1597,1577,1552,1524,1491,1455,1454,1490,1453,1413,1412,1368,1369,1414,1370,1415,1371,1372,1373,1417,1416,1457,1493,1456,1492,1525,1553,1598,1578,1554,1526,1527,1555,1599,1579,1459,1495,1494,1458,1418,1374,1375,1419,1327,1275,1159,1158,1219,1218,1157,1156,1217,1273,1274,1326,1325,1324,1272,1271,1323,1322,1270,1214,1213,1153,1154,1155,1215,1216,1091,1023,1022,1090,1089,1088,1020,1021,948,947,871,872,873,949,950,951,875,874,876,877,953,952,1024,1092,1093,1025,1026,1094,1095,1027,955,954,878,879,799,719,718,798,797,795,796,716,717,641,640,568,569,570,642,643,571,503,502,439,323,379,378,322,321,377,437,438,501,500,499,498,435,436,376,320,375,374,319,318,373,433,434,497,496,564,635,636,637,565,566,567,639,638,713,714,715,794,793,792,712,711,790,791,789,788,709,710,634,563,562,633,632,561,560,631,706,707,708,787,786,785,784,705,704,783,782,781,702,703,628,627,556,557,558,629,630,559,491,492,428,427,490,489,426,367,312,313,368,369,314,315,371,370,429,493,494,430,431,495,432,372,316,317,266,265,218,175,135,174,134,173,216,217,264,263,215,262,261,214,171,132,172,133,98,66,97,39,16,40,67,99,68,100,41,17,18,42,70,102,101,69,136,176,219,267,220,268,269,221,177,137,138,178,270,222,223,271,179,139,103,71,43,19,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631};
        // clang-format on
    }
    else if (grid_name() == "unstructured" && mpi::comm().size() == 1) {
        // clang-format off
        expected = {0,241,53,128,4,5,234,74,177,52,111,51,50,141,173,83,224,147,57,114,158,79,116,70,228,178,197,93,143,6,185,137,7,8,183,205,95,193,145,9,10,236,11,1,243,12,131,13,214,76,150,226,119,14,190,15,112,87,211,123,59,127,155,171,96,220,54,196,103,63,201,67,239,16,174,17,217,98,139,223,210,122,18,208,19,20,134,88,200,105,238,163,120,71,133,132,81,121,62,180,138,100,187,222,66,161,91,213,84,194,109,49,48,203,47,188,189,102,181,80,159,46,45,44,207,168,108,216,85,186,43,167,42,41,142,110,204,65,195,89,218,237,164,125,153,151,152,148,61,192,165,115,56,230,73,135,106,64,154,72,21,162,231,22,176,101,157,209,140,99,23,175,24,240,68,202,104,221,55,97,219,156,78,117,58,124,212,86,113,25,191,26,118,225,149,75,215,27,129,28,242,2,29,233,30,31,144,199,92,227,179,69,136,184,32,33,182,34,146,94,206,198,229,126,166,90,170,60,160,107,232,82,169,40,39,38,172,77,235,35,36,130,37,244,3};
        // clang-format on
    }
    else {
        Log::warning() << "No ordering will be tested" << std::endl;
    }
    test_reordering(reorder_config, expected);
}

CASE("test_reverse_cuthill_mckee_reordering") {
    auto reorder_config = option::type("reverse_cuthill_mckee");

    std::vector<idx_t> expected;
    if (grid_name() == "O16" && mpi::comm().size() == 1) {
        // clang-format off
        expected = {1579,1599,1555,1554,1578,1598,1527,1526,1525,1553,1577,1597,1495,1494,1493,1492,1524,1552,1576,1596,1459,1458,1457,1456,1455,1491,1523,1551,1575,1595,1419,1418,1417,1416,1415,1414,1454,1490,1522,1550,1574,1573,1594,1375,1374,1373,1372,1371,1370,1369,1413,1453,1489,1521,1549,1548,1547,1572,1593,1327,1326,1325,1324,1323,1322,1321,1320,1368,1412,1452,1488,1520,1519,1518,1517,1546,1571,1592,1275,1274,1273,1272,1271,1270,1269,1268,1267,1319,1367,1411,1451,1487,1486,1485,1484,1483,1516,1545,1570,1591,1219,1218,1217,1216,1215,1214,1213,1212,1211,1210,1266,1318,1366,1410,1450,1449,1448,1447,1446,1445,1482,1515,1544,1569,1590,1159,1095,1158,1157,1156,1155,1154,1153,1152,1151,1150,1149,1209,1148,1265,1317,1365,1409,1408,1407,1406,1405,1404,1589,1403,1444,1481,1514,1543,1568,1567,1027,955,799,879,1094,1026,1093,1092,1091,1090,1089,1088,1087,1086,1085,1084,1083,1208,1147,1264,1316,1364,1363,1362,1361,1360,1359,1358,1566,1588,1357,1402,1443,1480,1513,1542,1541,1540,878,954,719,798,1082,1025,953,1024,1023,1022,1021,1020,1019,1018,1017,1016,1015,1014,1013,1207,1146,1263,1315,1314,1313,1312,1311,1310,1309,1308,1539,1565,1587,1307,1356,1401,1442,1479,1512,1511,1510,1509,877,797,643,718,1081,1012,876,952,951,950,949,948,947,946,945,944,943,942,941,940,939,1206,1145,1262,1261,1260,1259,1258,1257,1256,1255,1254,1508,1538,1564,1586,1253,1306,1355,1400,1441,1478,1477,1476,1475,1474,796,717,571,642,1080,1011,938,874,875,872,873,870,871,868,869,866,867,864,865,862,863,861,1205,1144,1204,1203,1202,1201,1200,1199,1198,1197,1196,1473,1507,1537,1563,1585,1584,1583,1582,1580,1581,1195,1252,1305,1354,1399,1440,1439,1438,1437,1436,1435,795,716,641,503,570,1079,1010,860,937,794,793,792,791,790,789,788,787,786,785,784,783,782,780,781,1143,1078,1142,1141,1140,1139,1138,1137,1136,1135,1134,1133,1434,1472,1506,1536,1562,1561,1560,1559,1558,1556,1557,1194,1132,1251,1304,1353,1398,1397,1396,1395,1394,1393,1392,715,640,569,439,502,1009,936,859,779,714,713,712,711,710,709,708,707,706,705,704,703,702,701,1077,1008,1076,1075,1074,1073,1072,1071,1070,1069,1068,1067,1066,1391,1433,1471,1505,1535,1534,1533,1532,1531,1530,1528,1529,1193,1131,1250,1303,1352,1351,1350,1349,1348,1347,1346,1345,639,568,501,379,438,858,935,778,700,638,637,636,635,634,633,632,631,630,629,628,627,626,1065,1007,934,1006,1005,1004,1003,1002,1001,1000,999,998,997,996,995,1344,1390,1432,1470,1504,1503,1502,1501,1500,1499,1498,1496,1497,1192,1130,1249,1302,1301,1300,1299,1298,1297,1296,1295,1294,567,500,437,323,378,857,777,699,625,566,565,564,563,562,561,560,559,558,557,556,555,1064,994,856,933,932,931,930,929,928,927,926,925,924,923,922,921,920,1293,1343,1389,1431,1469,1468,1467,1466,1465,1464,1463,1462,1460,1461,1191,1129,1248,1247,1246,1245,1244,1243,1242,1241,1240,1239,499,436,377,271,322,776,698,624,554,498,497,496,495,494,493,492,491,490,489,488,1063,993,919,854,855,852,853,850,851,848,849,846,847,844,845,842,843,841,1238,1292,1342,1388,1430,1429,1428,1427,1426,1425,1424,1423,1422,1420,1421,1190,1128,1189,1188,1187,1186,1185,1184,1183,1182,1181,1180,435,376,321,223,270,775,697,623,553,487,434,433,432,431,430,429,428,427,426,425,1062,992,840,918,774,773,772,771,770,769,768,767,766,765,764,763,762,760,761,1179,1237,1291,1341,1387,1386,1385,1384,1383,1382,1381,1380,1379,1378,1376,1377,1127,1061,1126,1125,1124,1123,1122,1121,1120,1119,1118,1117,1116,375,320,269,179,222,696,622,552,486,424,374,373,372,371,370,369,368,367,366,991,917,839,759,695,694,693,692,691,690,689,688,687,686,685,684,683,682,1178,1236,1290,1340,1339,1338,1337,1336,1335,1334,1333,1332,1331,1330,1328,1329,1115,1060,990,1059,1058,1057,1056,1055,1054,1053,1052,1051,1050,1049,319,268,221,139,178,621,551,485,423,365,318,317,316,315,314,313,312,311,838,916,758,681,620,619,618,617,616,615,614,613,612,611,610,609,608,1177,1235,1289,1288,1287,1286,1285,1284,1283,1282,1281,1280,1279,1278,1276,1277,1048,1114,989,915,988,987,986,985,984,983,982,981,980,979,978,977,267,220,177,103,138,550,484,422,364,310,266,265,264,263,262,261,260,837,757,680,607,549,548,547,546,545,544,543,542,541,540,539,538,1176,1234,1233,1232,1231,1230,1229,1228,1227,1226,1225,1224,1223,1222,1220,1221,1047,976,1113,836,914,913,912,911,910,909,908,907,906,905,904,903,902,901,219,176,137,71,102,483,421,363,309,259,218,217,216,215,214,213,756,679,606,537,482,481,480,479,478,477,476,475,474,473,472,1175,1174,1173,1172,1171,1170,1169,1168,1167,1166,1165,1164,1163,1162,1160,1161,1046,975,1112,900,834,835,832,833,830,831,828,829,826,827,824,825,822,823,821,175,136,101,43,70,420,362,308,258,212,174,173,172,171,170,755,678,605,536,471,419,418,417,416,415,414,413,412,411,410,1111,1110,1109,1108,1107,1106,1105,1104,1103,1102,1101,1100,1099,1098,1096,1097,1045,974,1044,820,899,754,753,752,751,750,749,748,747,746,745,744,743,742,740,741,19,135,100,69,42,361,307,257,211,169,134,133,132,131,677,604,535,470,409,360,359,358,357,356,355,354,353,352,1043,1042,1041,1040,1039,1038,1037,1036,1035,1034,1033,1032,1031,1030,1028,1029,973,898,972,819,739,676,675,674,673,672,671,670,669,668,667,666,665,664,663,18,99,68,41,306,256,210,168,130,98,97,96,603,534,469,408,351,305,304,303,302,301,300,299,298,971,970,969,968,967,966,965,964,963,962,961,960,956,959,958,957,818,897,896,738,662,602,601,600,599,598,597,596,595,594,593,592,591,590,17,67,40,255,209,167,129,95,66,65,533,468,407,350,297,254,253,252,251,250,249,248,895,894,893,892,891,890,889,888,887,886,885,884,880,881,883,882,816,817,737,661,589,532,531,530,529,528,527,526,525,524,523,522,521,800,16,39,208,166,128,94,64,38,467,406,349,296,247,207,206,205,204,203,202,814,815,812,813,810,811,808,809,806,807,804,805,801,802,803,736,660,588,520,466,465,464,463,462,461,460,459,458,457,456,720,15,165,127,93,63,37,405,348,295,246,201,164,163,162,161,160,735,734,733,732,731,730,729,728,727,726,725,724,723,721,722,659,587,519,455,404,403,402,401,400,399,398,397,396,395,644,14,126,92,62,36,347,294,245,200,159,125,124,123,122,658,657,656,655,654,653,652,651,650,649,648,647,646,645,586,518,454,394,346,345,344,343,342,341,340,339,338,572,13,91,61,35,293,244,199,158,121,90,89,88,585,584,583,582,581,580,579,578,577,576,575,574,573,517,453,393,337,292,291,290,289,288,287,286,285,504,12,60,34,243,198,157,120,87,59,58,516,515,514,513,512,511,510,509,508,507,506,505,452,392,336,284,242,241,240,239,238,237,236,440,11,33,197,156,119,86,57,32,451,450,449,448,447,446,445,444,443,442,441,391,335,283,235,196,195,194,193,192,191,380,10,155,118,85,56,31,390,389,388,387,386,385,384,383,382,381,334,282,234,190,154,153,152,151,150,324,9,117,84,55,30,333,332,331,330,329,328,327,326,325,281,233,189,149,116,115,114,113,272,8,83,54,29,280,279,278,277,276,275,274,273,232,188,148,112,82,81,80,224,7,53,28,231,230,229,228,227,226,225,187,147,111,79,52,51,180,6,27,186,185,184,183,182,181,146,110,78,50,26,140,5,145,144,143,142,141,109,77,49,25,104,4,108,107,106,105,76,48,24,72,3,75,74,73,47,23,44,2,46,45,22,20,1,21,0,1631,1630,1629,1628,1627,1626,1625,1624,1623,1622,1621,1620,1619,1618,1617,1616,1615,1614,1613,1612,1611,1610,1609,1608,1607,1606,1605,1604,1603,1602,1601,1600};
        // clang-format on
    }
    else if (grid_name() == "unstructured" && mpi::comm().size() == 1) {
        // clang-format off
        expected = {2,242,29,28,33,32,31,30,184,129,233,34,35,27,215,182,136,144,227,92,75,3,36,244,94,146,235,26,118,229,69,199,179,149,37,130,198,206,77,160,25,86,191,24,117,126,58,225,38,172,60,107,232,212,113,22,175,23,166,156,78,124,39,170,169,82,21,231,68,240,176,99,90,153,97,219,202,40,164,142,204,237,20,19,162,101,140,209,125,55,104,41,65,110,218,89,208,18,134,72,157,64,151,221,135,17,42,195,167,148,152,88,217,122,200,154,106,73,230,16,174,43,85,186,216,165,61,210,98,223,105,56,115,238,81,15,239,112,44,108,168,192,207,63,139,163,133,120,132,121,14,190,13,12,87,67,211,45,80,159,201,103,71,196,180,62,181,1,243,119,214,131,150,226,123,46,102,96,220,54,138,91,213,189,11,76,236,193,59,155,47,188,171,187,222,161,84,10,95,145,205,127,79,48,203,66,100,109,9,183,158,70,116,49,194,114,141,8,137,228,57,178,50,83,224,173,7,93,185,147,197,51,111,6,143,74,52,177,5,234,128,53,4,241,0};
        // clang-format on
    }
    else {
        Log::warning() << "No ordering will be tested" << std::endl;
    }
    test_reordering(reorder_config, expected);
}

CASE("test_none_reordering") {
    auto reorder_config = option::type("none");
    test_reordering(reorder_config);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
