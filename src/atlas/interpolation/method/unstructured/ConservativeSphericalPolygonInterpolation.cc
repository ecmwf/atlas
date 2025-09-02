/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iomanip>
#include <vector>
#include <sys/stat.h> // for mkdir

#include "ConservativeSphericalPolygonInterpolation.h"

#include "eckit/log/ProgressTimer.h"

#include "atlas/grid.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/library/FloatingPointExceptions.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Topology.h"

#include "eckit/log/Bytes.h"

#define PRINT_BAD_POLYGONS 0

namespace atlas {
namespace interpolation {
namespace method {

MethodBuilder<ConservativeSphericalPolygonInterpolation> __builder("conservative-spherical-polygon");

using runtime::trace::StopWatch;
using util::ConvexSphericalPolygon;

namespace {

template <typename It>
std::string to_json(const It& begin, const It& end, int precision = 0) {
    std::stringstream ss;
    ss << "[\n";
    size_t size = std::distance(begin,end);
    size_t c=0;
    for( auto it = begin; it != end; ++it, ++c ) {
        ss << "  " << it->json(precision);
        if( c < size-1 ) {
            ss << ",\n";
        }
    }
    ss << "\n]";
    return ss.str();
}

template<typename ConvexSphericalPolygonContainer>
std::string to_json(const ConvexSphericalPolygonContainer& polygons, int precision = 0) {
    return to_json(polygons.begin(),polygons.end(),precision);
}

template <typename It>
std::string polygons_to_json(const It& begin, const It& end, int precision = 0) {
    std::stringstream ss;
    ss << "[\n";
    size_t size = std::distance(begin,end);
    size_t c=0;
    for( auto it = begin; it != end; ++it, ++c ) {
        ss << "  " << it->json(precision);
        if( c < size-1 ) {
            ss << ",\n";
        }
    }
    ss << "\n]";
    return ss.str();
}

template<typename ConvexSphericalPolygonContainer>
std::string polygons_to_json(const ConvexSphericalPolygonContainer& polygons, int precision = 0) {
    return polygons_to_json(polygons.begin(),polygons.end(),precision);
}

void dump_polygons_to_json( const ConvexSphericalPolygon& t_csp,
                    double pointsSameEPS,
                    const std::vector<std::tuple<ConvexSphericalPolygon, int>>& source_polygons,
                    const std::vector<idx_t>& source_polygons_considered_indices,
                    const std::string folder,
                    const std::string name) {
    std::vector<ConvexSphericalPolygon> csp_arr{ t_csp };
    std::vector<ConvexSphericalPolygon> csp_arr_intersecting {t_csp};
    std::vector<ConvexSphericalPolygon> intersections;
    int count = 1;
    for( auto& s_idx : source_polygons_considered_indices ) {
        auto s_csp = std::get<0>(source_polygons[s_idx]);
        csp_arr.emplace_back( s_csp );
        std::fstream file_plg_debug(folder + name + "_" + std::to_string(count++) + ".debug", std::ios::out);
        ConvexSphericalPolygon iplg = t_csp.intersect(s_csp, &file_plg_debug, pointsSameEPS);
        file_plg_debug.close();
        if (iplg) {
            if( iplg.area() > 0. ) {
                csp_arr_intersecting.emplace_back( s_csp );
                intersections.emplace_back( iplg );
            }
        }
    }
    double tot = 0.;
    Log::info().indent();
    std::fstream file_info(folder + name + ".info", std::ios::out);
    file_info << "List of intersection weights: " << std::endl;
    count = 1;
    for( auto& iplg : intersections ) {
        csp_arr.emplace_back( iplg );
        csp_arr_intersecting.emplace_back( iplg );
        tot += iplg.area() / t_csp.area();
        file_info << "\t" << count++ << ": " << iplg.area() / t_csp.area() << std::endl;
    }
    file_info << std::endl << name + ": " << 100. * tot << " % covered."<< std::endl << std::endl;
    file_info << "Target polygon + candidate source polygons + " << intersections.size() << " intersections in file:" << std::endl;
    file_info << "\t" << folder + name + ".candidates\n" << std::endl;
    std::fstream file_plg(folder + name + ".candidates", std::ios::out);
    file_plg << polygons_to_json(csp_arr, 16);
    file_plg.close();
    file_info << "Target polygon + " << csp_arr_intersecting.size() << " intersecting source polygon + " << intersections.size() << " intersections in file:" << std::endl;
    file_info << "\t" << folder + name + ".intersections\n" << std::endl;
    file_plg.open(folder + name + ".intersections", std::ios::out);
    file_plg << polygons_to_json(csp_arr_intersecting, 16);
    file_plg.close();
    Log::info().unindent();
}

constexpr double unit_sphere_area() {
    return 4. * M_PI;
}

template <typename T>
size_t memory_of(const std::vector<T>& vector) {
    return sizeof(T) * vector.capacity();
}

template <typename T>
size_t memory_of(const std::vector<std::vector<T>>& vector_of_vector) {
    size_t mem = 0;
    for (const auto& vector : vector_of_vector) {
        mem += memory_of(vector);
    }
    return mem;
}

size_t memory_of(
    const std::vector<ConservativeSphericalPolygonInterpolation::InterpolationParameters>& vector_of_params) {
    size_t mem = 0;
    for (const auto& params : vector_of_params) {
        mem += memory_of(params.cell_idx);
        mem += memory_of(params.centroids);
        mem += memory_of(params.weights);
        mem += memory_of(params.tgt_weights);
    }
    return mem;
}

Mesh extract_mesh(FunctionSpace fs) {
    if (functionspace::CellColumns(fs)) {
        return functionspace::CellColumns(fs).mesh();
    }
    else if (functionspace::NodeColumns(fs)) {
        return functionspace::NodeColumns(fs).mesh();
    }
    else {
        ATLAS_THROW_EXCEPTION("Cannot extract mesh from FunctionSpace" << fs.type());
    }
}

void sort_and_accumulate_triplets(std::vector<eckit::linalg::Triplet>& triplets) {
    ATLAS_TRACE();
    std::map<std::pair<int, int>, double> triplet_map;
    ATLAS_TRACE_SCOPE("accumulate in map")
    for (auto& triplet : triplets) {
        auto loc   = std::make_pair<int, int>(triplet.row(), triplet.col());
        auto entry = triplet_map.find(loc);
        if (entry == triplet_map.end()) {
            triplet_map[loc] = triplet.value();
        }
        else {
            entry->second += triplet.value();
        }
    }
    triplets.clear();
    ATLAS_TRACE_SCOPE("recontruct sorted vector from map")
    for (auto& triplet : triplet_map) {
        auto& row = triplet.first.first;
        auto& col = triplet.first.second;
        auto& val = triplet.second;
        triplets.emplace_back(row, col, val);
    }
}

int inside_vertices(const ConvexSphericalPolygon& plg1, const ConvexSphericalPolygon& plg2, int& pout) {
    int points_in = 0;
    pout          = 0;
    for (int j = 0; j < plg2.size(); j++) {
        int i = 0;
        for (; i < plg1.size(); i++) {
            int in   = (i != plg1.size() - 1) ? i + 1 : 0;
            auto gss = ConvexSphericalPolygon::GreatCircleSegment(plg1[i], plg1[in]);
            if (not gss.inLeftHemisphere(plg2[j], -std::numeric_limits<double>::epsilon())) {
                pout++;
                break;
            };
        }
        if (i == plg1.size()) {
            points_in++;
        }
    }
    ATLAS_ASSERT(points_in + pout == plg2.size());
    return points_in;
}

inline bool valid_point(idx_t node_idx, const array::ArrayView<int, 1>& node_flags) {
    return not util::Bitflags::view(node_flags(node_idx)).check(util::Topology::INVALID);
}

}  // namespace

ConservativeSphericalPolygonInterpolation::ConservativeSphericalPolygonInterpolation(const Config& config):
    Method(config) {
    config.get("validate", validate_ = false);
    config.get("order", order_ = 1);
    config.get("normalise_intersections", normalise_intersections_ = 0);
    config.get("matrix_free", matrix_free_ = false);
    config.get("src_cell_data", src_cell_data_ = true);
    config.get("tgt_cell_data", tgt_cell_data_ = true);


    config.get("statistics.timings", statistics_timings_ = false);
    config.get("statistics.intersection", statistics_intersection_ = false);
    config.get("statistics.conservation", statistics_conservation_ = false);
    if (statistics_conservation_ && ! statistics_timings_) {
        Log::info() << "statistics.conservation requested -> enabling statistics.timings";
        statistics_timings_ = true;
    }

    sharable_data_ = std::make_shared<Data>();
    cache_         = Cache(sharable_data_);
    data_          = sharable_data_.get();

    ATLAS_ASSERT(sharable_data_.use_count() == 2);
}

int ConservativeSphericalPolygonInterpolation::next_index(int current_index, int size, int offset) const {
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
    ATLAS_ASSERT(offset >= 0 && offset <= size);
    return (current_index < size - offset) ? current_index + offset : current_index + offset - size;
}
int ConservativeSphericalPolygonInterpolation::prev_index(int current_index, int size, int offset) const {
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
    ATLAS_ASSERT(offset >= 0 && offset <= size);
    return (current_index >= offset) ? current_index - offset : current_index - offset + size;
}

// get counter-clockwise sorted neighbours of a cell
std::vector<idx_t> ConservativeSphericalPolygonInterpolation::get_cell_neighbours(Mesh& mesh, idx_t cell) const {
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto& nodes_ll  = array::make_view<double, 2>(mesh.nodes().lonlat());
    const idx_t n_nodes   = cell2node.cols(cell);
    std::vector<idx_t> nbr_cells;
    nbr_cells.reserve(n_nodes);
    if (mesh.nodes().cell_connectivity().rows() == 0) {
        mesh::actions::build_node_to_cell_connectivity(mesh);
    }
    const auto& node2cell  = mesh.nodes().cell_connectivity();
    const auto n2c_missval = node2cell.missing_value();

    for (idx_t inode = 0; inode < n_nodes; ++inode) {
        idx_t node0             = cell2node(cell, inode);
        idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
        const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
        const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
        PointXYZ p0, p1;
        eckit::geometry::Sphere::convertSphericalToCartesian(1., p0_ll, p0);
        eckit::geometry::Sphere::convertSphericalToCartesian(1., p1_ll, p1);
        if (PointXYZ::norm(p0 - p1) < 1e-14) {
            continue;  // edge = point
        }
        bool still_search = true;  // still search the cell having vertices node0 & node1, not having index "cell"
        int n_cells0      = node2cell.cols(node0);
        int n_cells1      = node2cell.cols(node1);
        for (int icell0 = 0; still_search && icell0 < n_cells0; icell0++) {
            int cell0 = node2cell(node0, icell0);
            if (cell0 == cell) {
                continue;
            }
            for (int icell1 = 0; still_search && icell1 < n_cells1; icell1++) {
                int cell1 = node2cell(node1, icell1);
                if (cell0 == cell1 && cell0 != n2c_missval) {
                    nbr_cells.emplace_back(cell0);
                    still_search = false;
                }
            }
        }
    }
    return nbr_cells;
}

struct ConservativeSphericalPolygonInterpolation::Workspace {
    std::vector<idx_t> nbr_nodes_od;
    std::vector< std::array<idx_t,2> > cnodes;
};

// get cyclically sorted node neighbours without using edge connectivity
std::vector<idx_t> ConservativeSphericalPolygonInterpolation::get_node_neighbours(Mesh& mesh, idx_t node_id, Workspace& w) const {
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto node_flags = array::make_view<int, 1>(mesh.nodes().flags());
    if (mesh.nodes().cell_connectivity().rows() == 0) {
        mesh::actions::build_node_to_cell_connectivity(mesh);
    }
    const auto& node2cell = mesh.nodes().cell_connectivity();
    std::vector<idx_t> nbr_nodes;
    const int ncells = node2cell.cols(node_id);
    ATLAS_ASSERT(ncells > 0, "There is a node which does not connect to any cell");
    nbr_nodes.reserve(ncells + 1);
    w.cnodes.resize(ncells);
    w.nbr_nodes_od.clear();
    w.nbr_nodes_od.reserve(ncells + 1);
    for (idx_t icell = 0; icell < ncells; ++icell) {
        const idx_t cell = node2cell(node_id, icell);
        const int nnodes = cell2node.cols(cell);
        idx_t cnode      = 0;
        for (; cnode < nnodes; ++cnode) {
            if (node_id == cell2node(cell, cnode)) {
                break;
            }
        }
        w.cnodes[icell][0] = cell2node(cell, prev_index(cnode, nnodes));
        w.cnodes[icell][1] = cell2node(cell, next_index(cnode, nnodes));
    }
    if (ncells == 1) {
        nbr_nodes.emplace_back(w.cnodes[0][0]);
        nbr_nodes.emplace_back(w.cnodes[0][1]);
        return nbr_nodes;
    }
    // cycle one direction
    idx_t find = w.cnodes[0][1];
    idx_t prev = w.cnodes[0][0];
    nbr_nodes.emplace_back(prev);
    nbr_nodes.emplace_back(find);
    for (idx_t icycle = 0; nbr_nodes[0] != find;) {
        idx_t jcell = 0;
        for (; jcell < ncells; ++jcell) {
            idx_t ocell = (icycle + jcell + 1) % ncells;
            idx_t cand0 = w.cnodes[ocell][0];
            idx_t cand1 = w.cnodes[ocell][1];
            if (find == cand0 && prev != cand1) {
                if (cand1 == nbr_nodes[0]) {
                    return nbr_nodes;
                }
                nbr_nodes.emplace_back(cand1);
                prev = find;
                find = cand1;
                break;
            }
            if (find == cand1 && prev != cand0) {
                if (cand0 == nbr_nodes[0]) {
                    return nbr_nodes;
                }
                nbr_nodes.emplace_back(cand0);
                prev = find;
                find = cand0;
                break;
            }
        }
        if (jcell == ncells) {  // not found
            if (nbr_nodes[0] != find && find != nbr_nodes[nbr_nodes.size() - 1]) {
                nbr_nodes.emplace_back(find);
            }
            break;
        }
        else {
            icycle++;
        }
    }
    if (nbr_nodes[0] == find) {
        return nbr_nodes;
    }
    // cycle the oposite direction
    find = w.cnodes[0][0];
    prev = w.cnodes[0][1];
    w.nbr_nodes_od.emplace_back(prev);
    w.nbr_nodes_od.emplace_back(find);
    for (idx_t icycle = 0; w.nbr_nodes_od[0] != find;) {
        idx_t jcell = 0;
        for (; jcell < ncells; ++jcell) {
            idx_t ocell = (icycle + jcell + 1) % ncells;
            if (find == w.cnodes[ocell][0] && prev != w.cnodes[ocell][1]) {
                w.nbr_nodes_od.emplace_back(w.cnodes[ocell][1]);
                prev = find;
                find = w.cnodes[ocell][1];
                break;
            }
            if (find == w.cnodes[ocell][1] && prev != w.cnodes[ocell][0]) {
                w.nbr_nodes_od.emplace_back(w.cnodes[ocell][0]);
                prev = find;
                find = w.cnodes[ocell][0];
                break;
            }
        }
        if (jcell == ncells) {
            if (find != w.nbr_nodes_od[w.nbr_nodes_od.size() - 1]) {
                w.nbr_nodes_od.emplace_back(find);
            }
            break;
        }
        icycle++;
    }
    // put together
    int ow_size = w.nbr_nodes_od.size();
    for (int i = 0; i < ow_size - 2; i++) {
        if (valid_point(w.nbr_nodes_od[ow_size - 1 - i], node_flags)) {
            nbr_nodes.emplace_back(w.nbr_nodes_od[ow_size - 1 - i]);
        }
    }
    return nbr_nodes;
}

// Create polygons for cell-centred data. Here, the polygons are mesh cells
ConservativeSphericalPolygonInterpolation::CSPolygonArray
ConservativeSphericalPolygonInterpolation::get_polygons_celldata(FunctionSpace fs) const {
    CSPolygonArray cspolygons;
    auto mesh = extract_mesh(fs);
    const idx_t n_cells = mesh.cells().size();
    cspolygons.resize(n_cells);
    const auto& cell2node  = mesh.cells().node_connectivity();
    const auto lonlat      = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto cell_halo   = array::make_view<int, 1>(mesh.cells().halo());
    const auto& cell_flags = array::make_view<int, 1>(mesh.cells().flags());
    const auto& cell_part  = array::make_view<int, 1>(mesh.cells().partition());
    std::vector<PointLonLat> pts_ll;
    const int fs_halo = functionspace::CellColumns(fs).halo().size();
    for (idx_t cell = 0; cell < n_cells; ++cell) {
        if( cell_halo(cell) > fs_halo ) {
            continue;
        }
        int halo_type       = cell_halo(cell);
        const idx_t n_nodes = cell2node.cols(cell);
        pts_ll.clear();
        pts_ll.resize(n_nodes);
        for (idx_t jnode = 0; jnode < n_nodes; ++jnode) {
            idx_t inode   = cell2node(cell, jnode);
            pts_ll[jnode] = PointLonLat{lonlat(inode, 0), lonlat(inode, 1)};
        }
        const auto& bitflag = util::Bitflags::view(cell_flags(cell));
        if (bitflag.check(util::Topology::PERIODIC) and mpi::rank() == cell_part(cell)) {
            halo_type = -1;
        }
        std::get<0>(cspolygons[cell]) = ConvexSphericalPolygon(pts_ll);
        std::get<1>(cspolygons[cell]) = halo_type;
    }
    return cspolygons;
}

// Create polygons for cell-vertex data. Here, the polygons are subcells of mesh cells created as
// 	 (cell_centre, edge_centre, cell_vertex, edge_centre)
// additionally, subcell-to-node and node-to-subcells mapping are computed
ConservativeSphericalPolygonInterpolation::CSPolygonArray
ConservativeSphericalPolygonInterpolation::get_polygons_nodedata(FunctionSpace fs, std::vector<idx_t>& csp2node,
                                                                 std::vector<std::vector<idx_t>>& node2csp,
                                                                 std::array<double, 2>& errors) const {
    CSPolygonArray cspolygons;
    csp2node.clear();
    node2csp.clear();
    auto mesh = extract_mesh(fs);
    node2csp.resize(mesh.nodes().size());
    const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto cell_halo  = array::make_view<int, 1>(mesh.cells().halo());
    auto xyz2ll           = [](const atlas::PointXYZ& p_xyz) {
        PointLonLat p_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p_xyz, p_ll);
        return p_ll;
    };
    auto ll2xyz = [](const atlas::PointLonLat& p_ll) {
        PointXYZ p_xyz;
        eckit::geometry::Sphere::convertSphericalToCartesian(1., p_ll, p_xyz);
        return p_xyz;
    };
    const auto node_halo  = array::make_view<int, 1>(mesh.nodes().halo());
    const auto node_flags = array::make_view<int, 1>(mesh.nodes().flags());
    const auto node_part  = array::make_view<int, 1>(mesh.nodes().partition());

    idx_t cspol_id = 0;         // subpolygon enumeration
    errors         = {0., 0.};  // over/undershoots in creation of subpolygons
    const int fs_halo = functionspace::NodeColumns(fs).halo().size();

    std::vector<PointXYZ> pts_xyz;
    std::vector<PointLonLat> pts_ll;
    std::vector<int> pts_idx;

    for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
        if( cell_halo(cell) > fs_halo ) {
            continue;
        }
        ATLAS_ASSERT(cell < cell2node.rows());
        const idx_t n_nodes = cell2node.cols(cell);
        ATLAS_ASSERT(n_nodes > 2);
        PointXYZ cell_mid(0., 0., 0.);  // cell centre
        pts_xyz.clear();
        pts_ll.clear();
        pts_idx.clear();
        pts_xyz.reserve(n_nodes);
        pts_ll.reserve(n_nodes);
        pts_idx.reserve(n_nodes);
        for (idx_t inode = 0; inode < n_nodes; ++inode) {
            idx_t node0             = cell2node(cell, inode);
            idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
            const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
            const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
            PointXYZ p0             = ll2xyz(p0_ll);
            PointXYZ p1             = ll2xyz(p1_ll);
            if (PointXYZ::norm(p0 - p1) < 1e-14) {
                continue;  // skip this edge = a pole point
            }
            if (not valid_point(node0, node_flags)) {
                continue;
            }
            pts_xyz.emplace_back(p0);
            pts_ll.emplace_back(p0_ll);
            pts_idx.emplace_back(inode);
            cell_mid = cell_mid + p0;
            cell_mid = cell_mid + p1;
        }
        if (pts_xyz.size() < 3) {
            continue; // skip this cell
        }
        cell_mid                  = PointXYZ::div(cell_mid, PointXYZ::norm(cell_mid));
        PointLonLat cell_ll       = xyz2ll(cell_mid);
        double loc_csp_area_shoot = ConvexSphericalPolygon(pts_ll).area();
        // get ConvexSphericalPolygon for each valid edge
        int halo_type{0};
        for (int inode = 0; inode < pts_idx.size(); inode++) {
            int inode_n        = next_index(inode, pts_idx.size());
            idx_t node_n       = cell2node(cell, inode_n);
            PointXYZ iedge_mid = pts_xyz[inode] + pts_xyz[inode_n];
            iedge_mid          = PointXYZ::div(iedge_mid, PointXYZ::norm(iedge_mid));
            csp2node.emplace_back(node_n);
            node2csp[node_n].emplace_back(cspol_id);
            int inode_nn = next_index(inode_n, pts_idx.size());
            if (PointXYZ::norm(pts_xyz[inode_nn] - pts_xyz[inode_n]) < 1e-14) {
                ATLAS_THROW_EXCEPTION("Three cell vertices on a same great arc!");
            }
            PointXYZ jedge_mid;
            jedge_mid = pts_xyz[inode_nn] + pts_xyz[inode_n];
            jedge_mid = PointXYZ::div(jedge_mid, PointXYZ::norm(jedge_mid));
            std::array<PointLonLat,4> subpol_pts_ll;
            subpol_pts_ll[0] = cell_ll;
            subpol_pts_ll[1] = xyz2ll(iedge_mid);
            subpol_pts_ll[2] = pts_ll[inode_n];
            subpol_pts_ll[3] = xyz2ll(jedge_mid);
            halo_type    = node_halo(node_n);

            if (util::Bitflags::view(node_flags(node_n)).check(util::Topology::PERIODIC) and
                node_part(node_n) == mpi::rank()) {
                halo_type = -1;
            }

            ConvexSphericalPolygon cspi(subpol_pts_ll);
            loc_csp_area_shoot -= cspi.area();
            cspolygons.emplace_back(cspi, halo_type);
            cspol_id++;
        }
        if (halo_type == 0) {
            errors[0] += std::abs(loc_csp_area_shoot);
            errors[1] = std::max(std::abs(loc_csp_area_shoot), errors[1]);
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&errors[0], 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&errors[1], 1, eckit::mpi::max());
    }
    return cspolygons;
}

void ConservativeSphericalPolygonInterpolation::do_setup_impl(const Grid& src_grid, const Grid& tgt_grid) {
    ATLAS_TRACE("ConservativeMethod::do_setup( Grid, Grid )");
    ATLAS_ASSERT(src_grid);
    ATLAS_ASSERT(tgt_grid);

    tgt_fs_ = data_->tgt_fs_;
    src_fs_ = data_->src_fs_;

    if (not tgt_fs_) {
        auto tgt_mesh_config = tgt_grid.meshgenerator() | option::halo(0);
        ATLAS_TRACE_SCOPE("Generate target mesh") { tgt_mesh_ = MeshGenerator(tgt_mesh_config).generate(tgt_grid); }
        ATLAS_TRACE_SCOPE("Create target functionspace") {
            if (tgt_cell_data_) {
                tgt_fs_ = functionspace::CellColumns(tgt_mesh_, option::halo(0));
            }
            else {
                tgt_fs_ = functionspace::NodeColumns(tgt_mesh_, option::halo(1));
            }
        }
        sharable_data_->tgt_fs_ = tgt_fs_;
        ATLAS_ASSERT(data_->tgt_fs_);
    }

    if (not src_fs_) {
        auto src_mesh_config = src_grid.meshgenerator() | option::halo(2);
        ATLAS_TRACE_SCOPE("Generate source mesh") {
            if (mpi::size() > 1) {
                src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid, grid::MatchingPartitioner(tgt_mesh_));
            }
            else {
                src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid);
            }
        }
        ATLAS_TRACE_SCOPE("Create source functionspace") {
            if (src_cell_data_) {
                src_fs_ = functionspace::CellColumns(src_mesh_, option::halo(2));
            }
            else {
                src_fs_ = functionspace::NodeColumns(src_mesh_, option::halo(2));
            }
        }
        sharable_data_->src_fs_ = src_fs_;
        ATLAS_ASSERT(data_->tgt_fs_);
    }

    do_setup(src_fs_, tgt_fs_);
}


void ConservativeSphericalPolygonInterpolation::do_setup(const Grid& src_grid, const Grid& tgt_grid,
                                                         const interpolation::Cache& cache) {
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation::do_setup(Grid, Grid, Cache)");

    if (Cache(cache)) {
        Log::debug() << "Interpolation data found in cache -> no polygon intersections required" << std::endl;
        cache_ = Cache(cache);
        data_  = cache_.get();
        sharable_data_.reset();

        src_fs_ = data_->src_fs_;
        tgt_fs_ = data_->tgt_fs_;

        src_cell_data_ = functionspace::CellColumns(src_fs_);
        tgt_cell_data_ = functionspace::CellColumns(tgt_fs_);

        src_mesh_ = extract_mesh(src_fs_);
        tgt_mesh_ = extract_mesh(tgt_fs_);

        if (order_ == 1 && matrix_free_) {
            // We don't need to continue with setups required for first order matrix-free
            // such as mesh generation and functionspace creation.
            return;
        }
    }

    if (not matrix_free_) {
        auto matrix_cache = interpolation::MatrixCache(cache);
        if (matrix_cache) {
            if (matrix_cache.uid() == std::to_string(order_) || matrix_cache.uid().empty()) {
                Log::debug() << "Matrix found in cache -> no setup required at all" << std::endl;
                setMatrix(matrix_cache);
                return;
            }
        }
    }

    do_setup_impl(src_grid, tgt_grid);
}

void ConservativeSphericalPolygonInterpolation::do_setup(const FunctionSpace& source, const FunctionSpace& target,
                                                         const interpolation::Cache& cache) {
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation::do_setup(FunctionSpace, FunctionSpace, Cache)");

    if (not matrix_free_) {
        auto matrix_cache = interpolation::MatrixCache(cache);
        if (matrix_cache) {
            if (matrix_cache.uid() == std::to_string(order_) || matrix_cache.uid().empty()) {
                Log::debug() << "Matrix found in cache -> no setup required at all" << std::endl;
                setMatrix(matrix_cache);
                src_fs_ = source;
                tgt_fs_ = target;
                return;
            }
        }
        else {
            Log::warning() << "Could not find matrix in cache" << std::endl;
        }
    }

    if (Cache(cache)) {
        Log::debug() << "Interpolation data found in cache -> no polygon intersections required" << std::endl;
        cache_ = Cache(cache);
        data_  = cache_.get();
        sharable_data_.reset(new Data());
        tgt_fs_                 = target;
        src_fs_                 = source;
        src_cell_data_ = functionspace::CellColumns(src_fs_);
        tgt_cell_data_ = functionspace::CellColumns(tgt_fs_);

        src_mesh_ = extract_mesh(src_fs_);
        tgt_mesh_ = extract_mesh(tgt_fs_);

        if (order_ == 1 && matrix_free_) {
            // We don't need to continue with setups required for first order matrix-free
            // such as mesh generation and functionspace creation.
            return;
        }
    }

    do_setup(source, target);
}

void ConservativeSphericalPolygonInterpolation::do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs) {
    ATLAS_TRACE("ConservativeMethod::do_setup( FunctionSpace, FunctionSpace )");
    ATLAS_ASSERT(src_fs);
    ATLAS_ASSERT(tgt_fs);

    bool compute_cache = data_->src_points_.empty();

    if (not data_->tgt_fs_) {
        tgt_fs_                 = tgt_fs;
        sharable_data_->tgt_fs_ = tgt_fs_;
    }
    if (not data_->src_fs_) {
        src_fs_                 = src_fs;
        sharable_data_->src_fs_ = src_fs_;
    }

    src_cell_data_ = functionspace::CellColumns(src_fs_);
    tgt_cell_data_ = functionspace::CellColumns(tgt_fs_);

    src_mesh_ = extract_mesh(src_fs_);
    tgt_mesh_ = extract_mesh(tgt_fs_);

    {
        // we need src_halo_size >= 2, tgt_halo_size >= 0 for CellColumns
        // if target is NodeColumns, we need:
        //      tgt_halo_size >= 1 and
        //      src_halo_size large enough to cover the the target halo cells in the first row
        int halo_size = 0;
        src_mesh_.metadata().get("halo", halo_size);
        if (halo_size < 2) {
            Log::info() << "WARNING The halo size on source mesh should be at least 2.\n";
        }
        if (not tgt_cell_data_) {
            Log::info() << "WARNING The source cells should cover the first row of the target halos.\n";
        }
        halo_size = 0;
        tgt_mesh_.metadata().get("halo", halo_size);
        if (not tgt_cell_data_ and halo_size == 0) {
            Log::info() << "WARNING The halo size on target mesh should be at least 1 for the target NodeColumns.\n";
        }
    }
    CSPolygonArray src_csp;
    CSPolygonArray tgt_csp;
    if (compute_cache) {
        std::array<double, 2> errors = {0., 0.};
        ATLAS_TRACE("Get source polygons");
        StopWatch stopwatch;
        stopwatch.start();
        if (src_cell_data_) {
            src_csp = get_polygons_celldata(src_fs_);
        }
        else {
            src_csp =
                get_polygons_nodedata(src_fs_, sharable_data_->src_csp2node_, sharable_data_->src_node2csp_, errors);
        }
        stopwatch.stop();
        sharable_data_->timings.source_polygons_assembly = stopwatch.elapsed();
        remap_stat_.counts[Statistics::Counts::SRC_PLG]         = src_csp.size();
        remap_stat_.errors[Statistics::Errors::SRC_SUBPLG_L1]   = errors[0];
        remap_stat_.errors[Statistics::Errors::SRC_SUBPLG_LINF] = errors[1];

        errors = {0., 0.};
        ATLAS_TRACE("Get target polygons");
        stopwatch.start();
        if (tgt_cell_data_) {
            tgt_csp = get_polygons_celldata(tgt_fs_);
        }
        else {
            tgt_csp =
                get_polygons_nodedata(tgt_fs_, sharable_data_->tgt_csp2node_, sharable_data_->tgt_node2csp_, errors);
        }
        stopwatch.stop();
        sharable_data_->timings.target_polygons_assembly = stopwatch.elapsed();
        remap_stat_.counts[Statistics::Counts::TGT_PLG]         = tgt_csp.size();
        remap_stat_.errors[Statistics::Errors::TGT_SUBPLG_L1]   = errors[0];
        remap_stat_.errors[Statistics::Errors::TGT_SUBPLG_LINF] = errors[1];
    }
    else {
        remap_stat_.counts[Statistics::Counts::SRC_PLG]         = -1111;
        remap_stat_.errors[Statistics::Errors::SRC_SUBPLG_L1]   = -1111.;
        remap_stat_.errors[Statistics::Errors::SRC_SUBPLG_LINF] = -1111.;
        remap_stat_.counts[Statistics::Counts::TGT_PLG]         = -1111;
        remap_stat_.errors[Statistics::Errors::TGT_SUBPLG_L1]   = -1111.;
        remap_stat_.errors[Statistics::Errors::TGT_SUBPLG_LINF] = -1111.;
    }

    n_spoints_ = src_fs_.size();
    n_tpoints_ = tgt_fs_.size();

    if (compute_cache) {
        intersect_polygons(src_csp, tgt_csp);

        auto& src_points_ = sharable_data_->src_points_;
        src_points_.resize(n_spoints_);
        sharable_data_->src_areas_.resize(n_spoints_);
        auto& src_areas_v = sharable_data_->src_areas_;
        {
            ATLAS_TRACE("Store src_areas and src_point");
            if (src_cell_data_) {
                for (idx_t spt = 0; spt < n_spoints_; ++spt) {
                    const auto& s_csp = std::get<0>(src_csp[spt]);
                    src_points_[spt]  = s_csp.centroid();
                    src_areas_v[spt]  = s_csp.area();
                }
            }
            else {
                auto& src_node2csp_ = sharable_data_->src_node2csp_;
                const auto lonlat   = array::make_view<double, 2>(src_mesh_.nodes().lonlat());
                for (idx_t spt = 0; spt < n_spoints_; ++spt) {
                    if (src_node2csp_[spt].size() == 0) {
                        // this is a node to which no subpolygon is associated
                        // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                        auto p = PointLonLat{lonlat(spt, 0), lonlat(spt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, src_points_[spt]);
                    }
                    else {
                        // .. in the other case, start computing the barycentre
                        src_points_[spt] = PointXYZ{0., 0., 0.};
                    }
                    src_areas_v[spt] = 0.;
                    for (idx_t isubcell = 0; isubcell < src_node2csp_[spt].size(); ++isubcell) {
                        idx_t subcell     = src_node2csp_[spt][isubcell];
                        const auto& s_csp = std::get<0>(src_csp[subcell]);
                        src_areas_v[spt] += s_csp.area();
                        src_points_[spt] = src_points_[spt] + PointXYZ::mul(s_csp.centroid(), s_csp.area());
                    }
                    double src_point_norm = PointXYZ::norm(src_points_[spt]);
                    if (src_point_norm == 0.) {
                        ATLAS_DEBUG_VAR(src_point_norm);
                        ATLAS_DEBUG_VAR(src_points_[spt]);
                        ATLAS_DEBUG_VAR(src_node2csp_[spt].size());
                        for (idx_t isubcell = 0; isubcell < src_node2csp_[spt].size(); ++isubcell) {
                            idx_t subcell     = src_node2csp_[spt][isubcell];
                            ATLAS_DEBUG_VAR(subcell);
                            const auto& s_csp = std::get<0>(src_csp[subcell]);
                            s_csp.print(Log::info());
                            Log::info() << std::endl;
                            src_areas_v[spt] += s_csp.area();
                            ATLAS_DEBUG_VAR(s_csp.area());
                            ATLAS_DEBUG_VAR(s_csp.centroid());
                            src_points_[spt] = src_points_[spt] + PointXYZ::mul(s_csp.centroid(), s_csp.area());
                        }
                        Log::info().flush();
                        // something went wrong, improvise
                        src_point_norm = 1.;
                    }
                    src_points_[spt] = PointXYZ::div(src_points_[spt], src_point_norm);
                }
            }
        }
        auto& tgt_points_ = sharable_data_->tgt_points_;
        tgt_points_.resize(n_tpoints_);
        sharable_data_->tgt_areas_.resize(n_tpoints_);
        auto& tgt_areas_v = sharable_data_->tgt_areas_;
        {
            ATLAS_TRACE("Store src_areas and src_point");
            if (tgt_cell_data_) {
                for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
                    const auto& t_csp = std::get<0>(tgt_csp[tpt]);
                    tgt_points_[tpt]  = t_csp.centroid();
                    tgt_areas_v[tpt]  = t_csp.area();
                }
            }
            else {
                auto& tgt_node2csp_ = sharable_data_->tgt_node2csp_;
                const auto lonlat   = array::make_view<double, 2>(tgt_mesh_.nodes().lonlat());
                for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
                    if (tgt_node2csp_[tpt].size() == 0) {
                        // this is a node to which no subpolygon is associated
                        // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                        auto p = PointLonLat{lonlat(tpt, 0), lonlat(tpt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, tgt_points_[tpt]);
                    }
                    else {
                        // .. in the other case, start computing the barycentre
                        tgt_points_[tpt] = PointXYZ{0., 0., 0.};
                    }
                    tgt_areas_v[tpt] = 0.;
                    for (idx_t isubcell = 0; isubcell < tgt_node2csp_[tpt].size(); ++isubcell) {
                        idx_t subcell     = tgt_node2csp_[tpt][isubcell];
                        const auto& t_csp = std::get<0>(tgt_csp[subcell]);
                        tgt_areas_v[tpt] += t_csp.area();
                        tgt_points_[tpt] = tgt_points_[tpt] + PointXYZ::mul(t_csp.centroid(), t_csp.area());
                    }
                    double tgt_point_norm = PointXYZ::norm(tgt_points_[tpt]);
                    if (tgt_point_norm == 0.) {
                        for (idx_t isubcell = 0; isubcell < tgt_node2csp_[tpt].size(); ++isubcell) {
                            idx_t subcell     = tgt_node2csp_[tpt][isubcell];
                            ATLAS_DEBUG_VAR(subcell);
                            const auto& t_csp = std::get<0>(tgt_csp[subcell]);
                            t_csp.print(Log::info());
                            Log::info() << std::endl;
                            tgt_areas_v[tpt] += t_csp.area();
                            ATLAS_DEBUG_VAR(t_csp.area());
                            ATLAS_DEBUG_VAR(t_csp.centroid());
                            tgt_points_[tpt] = tgt_points_[tpt] + PointXYZ::mul(t_csp.centroid(), t_csp.area());
                        }
                        Log::info().flush();
                        // something went wrong, improvise
                        tgt_point_norm = 1.;
                    }
                    tgt_points_[tpt] = PointXYZ::div(tgt_points_[tpt], tgt_point_norm);
                }
            }
        }
    }

    if (not matrix_free_) {
        StopWatch stopwatch;
        stopwatch.start();
        switch (order_) {
            case 1: {
                setMatrix(n_tpoints_, n_spoints_, compute_1st_order_triplets(), "1");
                break;
            }
            case 2: {
                setMatrix(n_tpoints_, n_spoints_, compute_2nd_order_triplets(), "2");
                break;
            }
            default: {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        stopwatch.stop();
        if (compute_cache) {
            sharable_data_->timings.matrix_assembly = stopwatch.elapsed();
        }
    }

    data_->print(Log::debug());

    if (statistics_intersection_) {
        setup_stat();
    }
}

void ConservativeSphericalPolygonInterpolation::intersect_polygons(const CSPolygonArray& src_csp,
                                                                   const CSPolygonArray& tgt_csp) {
    ATLAS_TRACE();
    auto& timings = sharable_data_->timings;
    StopWatch stopwatch;
    stopwatch.start();
    util::KDTree<idx_t> kdt_search;
    double max_tgtcell_rad = 0.;
    const int tgt_halo_intersection_depth = (tgt_cell_data_ ? 0 : 1); // if target NodeColumns, one target halo required for subcells around target nodes
    ATLAS_TRACE_SCOPE("build kd-tree for target polygons") {
        kdt_search.reserve(tgt_csp.size());
        for (idx_t jcell = 0; jcell < tgt_csp.size(); ++jcell) {
            auto tgt_halo_type = std::get<1>(tgt_csp[jcell]);
            if (tgt_halo_type <= tgt_halo_intersection_depth) { // and tgt_halo_type != -1) {
                const auto& t_csp = std::get<0>(tgt_csp[jcell]);
                kdt_search.insert(t_csp.centroid(), jcell);
                max_tgtcell_rad = std::max(max_tgtcell_rad, t_csp.radius());
            }
        }
        kdt_search.build();
    }
    stopwatch.stop();
    timings.target_kdtree_assembly = stopwatch.elapsed();

    StopWatch stopwatch_src_already_in;
    StopWatch stopwatch_kdtree_search;
    StopWatch stopwatch_polygon_intersections;

    stopwatch_src_already_in.start();


    // needed for intersect_polygons only, merely for detecting duplicate points
    // Treshold at which points are considered same
    double compare_pointxyz_eps = 1.e8 * std::numeric_limits<double>::epsilon();
    const char* ATLAS_COMPAREPOINTXYZ_EPS_FACTOR = ::getenv("ATLAS_COMPAREPOINTXYZ_EPS_FACTOR");
    if (ATLAS_COMPAREPOINTXYZ_EPS_FACTOR != nullptr) {
        compare_pointxyz_eps = std::atof(ATLAS_COMPAREPOINTXYZ_EPS_FACTOR) * std::numeric_limits<double>::epsilon();
    }

    auto compare_pointxyz = [eps=compare_pointxyz_eps] (const PointXYZ& f, const PointXYZ& s) -> bool {
        if (f[0] < s[0] - eps) {
            return true;
        }
        else if (std::abs(f[0] - s[0]) < eps) {
            if (f[1] < s[1] - eps) {
                return true;
            }
            else if (std::abs(f[1] - s[1]) < eps) {
                if (f[2] < s[2] - eps) {
                    return true;
                }
            }
        }
        return false;
    };

    std::set<PointXYZ, decltype(compare_pointxyz)> src_cent(compare_pointxyz);
    auto src_already_in = [&](const PointXYZ& halo_cent) {
        if (src_cent.find(halo_cent) == src_cent.end()) {
            src_cent.insert(halo_cent);
            return false;
        }
        return true;
    };
    stopwatch_src_already_in.stop();

    enum MeshSizeId
    {
        SRC,
        TGT,
        SRC_TGT_INTERSECT,
        SRC_NONINTERSECT
    };
    std::array<size_t, 4> num_pol{0, 0, 0, 0};
    enum AreaCoverageId
    {
        TOTAL_SRC,
        MAX_SRC
    };
    std::array<double, 2> area_coverage{0., 0.};
    auto& src_iparam_ = sharable_data_->src_iparam_;
    src_iparam_.resize(src_csp.size());

    std::vector<InterpolationParameters> tgt_iparam;  // only used for debugging
    if (validate_) {
        tgt_iparam.resize(tgt_csp.size());
    }

    // the worst target polygon coverage for analysis of intersection
    std::pair<idx_t, double> worst_tgt_overcover;
    std::pair<idx_t, double> worst_tgt_undercover;
    worst_tgt_overcover.first = -1;
    worst_tgt_overcover.second = -1.;
    worst_tgt_undercover.first = -1;
    worst_tgt_undercover.second = M_PI;

    // NOTE: polygon vertex points at distance < pointsSameEPS will be replaced with one point 
    constexpr double pointsSameEPS = 5.e6 * std::numeric_limits<double>::epsilon();

    bool fpe_for_polygon = ConvexSphericalPolygon::fpe();
    ConvexSphericalPolygon::fpe(false);
    bool fpe_disabled = fpe_for_polygon ? atlas::library::disable_floating_point_exception(FE_INVALID) : false;
    auto restore_fpe = [fpe_disabled, fpe_for_polygon] {
        if (fpe_disabled) {
            atlas::library::enable_floating_point_exception(FE_INVALID);
        }
        ConvexSphericalPolygon::fpe(fpe_for_polygon);
    };

    eckit::Channel blackhole;
    Log::debug() << "Find intersections for " << src_csp.size() << " source polygons" << std::endl;

    {

    eckit::ProgressTimer progress("Intersecting polygons ", 100, " percent", double(10),
                                  src_csp.size() > 50 ? Log::info() : blackhole);
    float last_progress_percent = 0.00;
    for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
        stopwatch_src_already_in.start();
        bool already_in = src_already_in((std::get<0>(src_csp[scell]).centroid()));
        stopwatch_src_already_in.stop();

        if (not already_in) {
            const auto& s_csp       = std::get<0>(src_csp[scell]);
            const double s_csp_area = s_csp.area();
            if (s_csp_area == 0.) {
                Log::warning() << "Skipping source polygon " << scell << " with area = 0" << std::endl;
                continue;
            }
            double src_cover_area   = 0.;
            if (statistics_timings_) { stopwatch_kdtree_search.start(); }
            auto tgt_cells = kdt_search.closestPointsWithinRadius(s_csp.centroid(), s_csp.radius() + max_tgtcell_rad);
            if (statistics_timings_) { stopwatch_kdtree_search.stop(); }
            for (idx_t ttcell = 0; ttcell < tgt_cells.size(); ++ttcell) {
                auto tcell        = tgt_cells[ttcell].payload();
                const auto& t_csp = std::get<0>(tgt_csp[tcell]);
                if (statistics_timings_) { stopwatch_polygon_intersections.start(); }
                ConvexSphericalPolygon csp_i = s_csp.intersect(t_csp, nullptr, pointsSameEPS);
                double csp_i_area            = csp_i.area();
                if (statistics_timings_) { stopwatch_polygon_intersections.stop(); }
                if (validate_) {
                    int pout;
                    // TODO: this can be removed soon
                    if (inside_vertices(s_csp, t_csp, pout) > 2 && csp_i.area() < 3e-16) {
                        dump_intersection("Zero area intersections with inside_vertices", s_csp, tgt_csp, tgt_cells);
                    }
                    // TODO: assuming intersector search works fine, this should be move under "if (csp_i_area > 0)"
                    tgt_iparam[tcell].cell_idx.emplace_back(scell);
                    tgt_iparam[tcell].tgt_weights.emplace_back(csp_i_area);
                }
                if (csp_i_area > 0) {
                    src_iparam_[scell].cell_idx.emplace_back(tcell);
                    src_iparam_[scell].weights.emplace_back(csp_i_area);
                    double target_weight = csp_i_area / t_csp.area();
                    src_iparam_[scell].tgt_weights.emplace_back(target_weight); // TODO: tgt_weights vector should be removed for the sake of highres
                    if (order_ == 2 or not matrix_free_ or not matrixAllocated()) {
                        src_iparam_[scell].centroids.emplace_back(csp_i.centroid());
                    }
                    src_cover_area += csp_i_area;

                    if (validate_) {
                        // this check is concerned with accuracy of polygon intersections
                        if (target_weight > 1.1) {
                            dump_intersection("Intersection larger than target", s_csp, tgt_csp, tgt_cells);
                        }
                        if (csp_i_area / s_csp_area > 1.1) {
                            dump_intersection("Intersection larger than source", s_csp, tgt_csp, tgt_cells);
                        }
                    }
                }
            }
            const double src_cover_err         = std::abs(s_csp_area - src_cover_area);
            const double src_cover_err_percent = 100. * src_cover_err / s_csp_area;
            if (src_cover_err_percent > 0.1 and std::get<1>(src_csp[scell]) == 0) {
                // NOTE: source cell at process boundary will not be covered by target cells, skip them
                // TODO: mark these source cells beforehand and compute error in them among the processes
                if (validate_ and mpi::size() == 1) {
                    dump_intersection("Source cell not exactly covered", s_csp, tgt_csp, tgt_cells);
                    if (statistics_intersection_) {
                        area_coverage[TOTAL_SRC] += src_cover_err;
                        area_coverage[MAX_SRC] = std::max(area_coverage[MAX_SRC], src_cover_err);
                    }
                }
            }
            if (src_iparam_[scell].cell_idx.size() == 0 and statistics_intersection_) {
                num_pol[SRC_NONINTERSECT]++;
            }
            if (normalise_intersections_ && src_cover_err_percent < 1.) {
                double wfactor = s_csp.area() / (src_cover_area > 0. ? src_cover_area : 1.);
                for (idx_t i = 0; i < src_iparam_[scell].weights.size(); i++) {
                    src_iparam_[scell].weights[i] *= wfactor;
                    src_iparam_[scell].tgt_weights[i] *= wfactor;
                }
            }
            if (statistics_intersection_) {
                num_pol[SRC_TGT_INTERSECT] += src_iparam_[scell].weights.size();
            }
        } // already in
        if ( double(scell) / double(src_csp.size()) > last_progress_percent ) {
            last_progress_percent += 0.01;
            ++progress;
        }
    }

    }

    restore_fpe();

    timings.polygon_intersections  = stopwatch_polygon_intersections.elapsed();
    timings.target_kdtree_search   = stopwatch_kdtree_search.elapsed();
    timings.source_polygons_filter = stopwatch_src_already_in.elapsed();
    num_pol[SRC]                   = src_csp.size();
    num_pol[TGT]                   = tgt_csp.size();
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(num_pol.data(), num_pol.size(), eckit::mpi::sum());
    }
    remap_stat_.counts[Statistics::Counts::INT_PLG]   = num_pol[SRC_TGT_INTERSECT];
    remap_stat_.counts[Statistics::Counts::UNCVR_SRC] = num_pol[SRC_NONINTERSECT];

    const std::string polygon_intersection_folder = "polygon_intersection/";
    if (validate_ && mpi::rank() == 0) {
        if (mkdir(polygon_intersection_folder.c_str(), 0777) != 0) {
            Log::info() << "WARNING Polygon intersection relevant information in is the folder \e[1mpolygon_intersection\e[0m." << std::endl;
        }
        else {
            Log::info() << "WARNING Could not create the folder \e[1mpolygon_intersection\e[0m." << std::endl;
        }
    }

    if (validate_) {
        double geo_err_l1   = 0.;
        double geo_err_linf = -1.;
        for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
            if (std::get<1>(src_csp[scell]) != 0 ) {
                // skip periodic & halo cells
                continue;
            }
            double diff_cell = std::get<0>(src_csp[scell]).area();
            for (idx_t icell = 0; icell < src_iparam_[scell].weights.size(); ++icell) {
                diff_cell -= src_iparam_[scell].weights[icell];
            }
            geo_err_l1 += std::abs(diff_cell);
            geo_err_linf = std::max(geo_err_linf, std::abs(diff_cell));
        }
        ATLAS_TRACE_MPI(ALLREDUCE) {
            mpi::comm().allReduceInPlace(geo_err_l1, eckit::mpi::sum());
            mpi::comm().allReduceInPlace(geo_err_linf, eckit::mpi::max());
        }
        if (mpi::size() == 1) {
            remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_L1]   = geo_err_l1 / unit_sphere_area();
            remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_LINF] = geo_err_linf;
        }
        else {
            remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_L1]    = -1111.;
            remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_LINF]  = -1111.;
        }

        std::fstream polygon_intersection_info("polygon_intersection", std::ios::out); 

        geo_err_l1   = 0.;
        geo_err_linf = -1.;
        for (idx_t tcell = 0; tcell < tgt_csp.size(); ++tcell) {
            if (std::get<1>(tgt_csp[tcell]) != 0) {
                // skip periodic & halo cells
                continue;
            }
            const auto& t_csp     = std::get<0>(tgt_csp[tcell]);


            double tgt_cover_area = 0.;
            const auto& tiparam   = tgt_iparam[tcell];

#if PRINT_BAD_POLYGONS
            // dump polygons in json format
            idx_t tcell_printout = 120;
            if (tcell == tcell_printout) {
                dump_polygons_to_json(t_csp, 1.e-14, src_csp, tiparam.cell_idx, "polygon_dump", "tcell" + std::to_string(tcell_printout));
            }
#endif
            for (idx_t icell = 0; icell < tiparam.cell_idx.size(); ++icell) {
                tgt_cover_area += tiparam.tgt_weights[icell];
            }
            /*
            // TODO: normalise to target cell
            double normm = t_csp.area() / (tgt_cover_area > 0. ? tgt_cover_area : t_csp.area());
            for (idx_t icell = 0; icell < tiparam.cell_idx.size(); ++icell) {
                idx_t scell = tiparam.cell_idx[icell];
                auto siparam = src_iparam_[scell];
                size_t tgt_intersectors = siparam.cell_idx.size();
                for (idx_t sicell = 0; sicell < tgt_intersectors; sicell++ ) {
                    if (siparam.cell_idx[icell] == tcell) {;
                        siparam.weights[icell] *= normm;
                        siparam.tgt_weights[icell] *= normm;
                    }
                }
            }
            */
            double diff_cell = tgt_cover_area - t_csp.area();
            geo_err_l1 += std::abs(diff_cell);
            geo_err_linf = std::max(geo_err_linf, std::abs(diff_cell));
            const double tgt_cover_err = 100. * diff_cell / t_csp.area();
            if (worst_tgt_overcover.second < tgt_cover_err) {
                worst_tgt_overcover.second = tgt_cover_err;;
                worst_tgt_overcover.first = tcell;
            }
            if (worst_tgt_undercover.second > tgt_cover_err) {
                worst_tgt_undercover.second = tgt_cover_err;;
                worst_tgt_undercover.first = tcell;
            }
            if (std::abs(tgt_cover_err) > 0.5) {
                PointLonLat centre_ll;
                eckit::geometry::Sphere::convertCartesianToSpherical(1., t_csp.centroid(), centre_ll);
                polygon_intersection_info << "WARNING tgt cell " << tcell << " over-covering: \e[1m" << tgt_cover_err << "\e[0m %, cell-centre: "
                    << centre_ll <<"\n";
                polygon_intersection_info << "source indices: " << tiparam.cell_idx << std::endl;
                dump_intersection("Target cell not exaclty covered", t_csp, src_csp, tiparam.cell_idx);
                //dump_polygons_to_json(t_csp, src_csp, tiparam.cell_idx, "bad_polygon", 1.e-16);
            }
        }
        ATLAS_TRACE_MPI(ALLREDUCE) {
            mpi::comm().allReduceInPlace(geo_err_l1, eckit::mpi::sum());
            mpi::comm().allReduceInPlace(geo_err_linf, eckit::mpi::max());
        }
        remap_stat_.errors[Statistics::Errors::TGT_INTERSECTPLG_L1]   = geo_err_l1 / unit_sphere_area();
        remap_stat_.errors[Statistics::Errors::TGT_INTERSECTPLG_LINF] = geo_err_linf;
    }
    else {
        remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_L1]    = -1111.;
        remap_stat_.errors[Statistics::Errors::SRC_INTERSECTPLG_LINF]  = -1111.;
        remap_stat_.errors[Statistics::Errors::TGT_INTERSECTPLG_L1]   = -1111.;
        remap_stat_.errors[Statistics::Errors::TGT_INTERSECTPLG_LINF] = -1111.;
    }

    if (validate_) {
        std::vector<idx_t> first(mpi::comm().size());
        std::vector<double> second(mpi::comm().size());
        ATLAS_TRACE_MPI(ALLGATHER) {
            mpi::comm().allGather(worst_tgt_overcover.first, first.begin(), first.end());
            mpi::comm().allGather(worst_tgt_overcover.second, second.begin(), second.end());
        }
        auto max_over = std::max_element(second.begin(), second.end());
        auto rank_over = std::distance(second.begin(), max_over);
        Log::info() << "WARNING The worst target polygon over-coveraging: \e[1m" 
            << *max_over
            << "\e[0m %. For details, check the file: worst_target_cell_overcover.info " << std::endl;
        if (rank_over == mpi::rank()) {
            auto tcell = worst_tgt_overcover.first;
            dump_polygons_to_json(std::get<0>(tgt_csp[tcell]), pointsSameEPS, src_csp, tgt_iparam[tcell].cell_idx, polygon_intersection_folder, "worst_target_cell_overcover");
        }
        ATLAS_TRACE_MPI(ALLGATHER) {
            mpi::comm().allGather(worst_tgt_undercover.first, first.begin(), first.end());
            mpi::comm().allGather(worst_tgt_undercover.second, second.begin(), second.end());
        }
        auto min_under = std::min_element(second.begin(), second.end());
        auto rank_under = std::distance(second.begin(), min_under);
        Log::info() << "WARNING The worst target polygon under-coveraging: \e[1m" 
            << *min_under
            << "\e[0m %. For details, check the file: worst_target_cell_undercover.info " << std::endl;

        if (rank_under == mpi::rank()) {
            auto tcell = worst_tgt_undercover.first;
            dump_polygons_to_json(std::get<0>(tgt_csp[tcell]), pointsSameEPS, src_csp, tgt_iparam[tcell].cell_idx, polygon_intersection_folder, "worst_target_cell_undercover");
        }
    }
}

ConservativeSphericalPolygonInterpolation::Triplets ConservativeSphericalPolygonInterpolation::compute_1st_order_triplets() {
    ATLAS_TRACE("ConservativeMethod::setup: build cons-1 interpolant matrix");
    ATLAS_ASSERT(not matrix_free_);
    Triplets triplets;
    size_t triplets_size    = 0;
    const auto& src_iparam_ = data_->src_iparam_;
    // determine the size of array of triplets used to define the sparse matrix
    if (src_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            triplets_size += src_iparam_[scell].cell_idx.size();
        }
    }
    else {
        auto& src_node2csp_ = data_->src_node2csp_;
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += src_iparam_[subcell].tgt_weights.size();
            }
        }
    }
    triplets.reserve(triplets_size);
    // assemble triplets to define the sparse matrix
    const auto& tgt_areas_v = data_->tgt_areas_;
    if (src_cell_data_ && tgt_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = src_iparam_[scell];
            for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                idx_t tcell = iparam.cell_idx[icell];
                triplets.emplace_back(tcell, scell, iparam.tgt_weights[icell]);
            }
        }
    }
    else if (not src_cell_data_ && tgt_cell_data_) {
        auto& src_node2csp_ = data_->src_node2csp_;
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = src_iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    ATLAS_ASSERT(tcell < n_tpoints_);
                    triplets.emplace_back(tcell, snode, iparam.tgt_weights[icell]);
                }
            }
        }
    }
    else if (src_cell_data_ && not tgt_cell_data_) {
        auto& tgt_csp2node_ = data_->tgt_csp2node_;
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = src_iparam_[scell];
            for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                idx_t tcell            = iparam.cell_idx[icell];
                idx_t tnode            = tgt_csp2node_[tcell];
                if (tnode >= n_tpoints_) {
                    Log::info() << "tnode, n_tpoints = " << tnode << ", " << n_tpoints_ << std::endl;
                    ATLAS_ASSERT(false);
                }
                double inv_node_weight = (tgt_areas_v[tnode] > 0. ? 1. / tgt_areas_v[tnode] : 0.);
                triplets.emplace_back(tnode, scell, iparam.weights[icell] * inv_node_weight);
            }
        }
    }
    else if (not src_cell_data_ && not tgt_cell_data_) {
        auto& src_node2csp_ = data_->src_node2csp_;
        auto& tgt_csp2node_ = data_->tgt_csp2node_;
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = src_iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    idx_t tcell            = iparam.cell_idx[icell];
                    idx_t tnode            = tgt_csp2node_[tcell];
                    ATLAS_ASSERT(tnode < n_tpoints_);
                    double inv_node_weight = (tgt_areas_v[tnode] > 0. ? 1. / tgt_areas_v[tnode] : 0.);
                    triplets.emplace_back(tnode, snode, iparam.weights[icell] * inv_node_weight);
                    if( tnode >= n_tpoints_) {
                        Log::info() << tnode << " = tnode, " << n_tpoints_ << " = n_tpoints\n";;
                        ATLAS_ASSERT(false);
                    }
                }
            }
        }
    }
    sort_and_accumulate_triplets(triplets);  // Very expensive!!! (90% of this routine). We need to avoid it

//     if (validate_) {
//         std::vector<double> weight_sum(n_tpoints_);
//         for( auto& triplet : triplets ) {
//             weight_sum[triplet.row()] += triplet.value();
//         }
//         if (order_ == 1) {
//             // first order should not give overshoots
//             double eps = 1e4 * std::numeric_limits<double>::epsilon();
//             for( auto& triplet : triplets ) {
//                 if (triplet.value() > 1. + eps or triplet.value() < -eps) {
//                     Log::info() << "target point " << triplet.row() << " weight: " << triplet.value() << std::endl;
//                 }
//             }
//         }
//         for( size_t row=0; row < n_tpoints_; ++row ) {
//             if (std::abs(weight_sum[row] - 1.) > 1e-11) {
//                 Log::info() << "target weight in row " << row << " differs from 1 by " << std::abs(weight_sum[row] - 1.) << std::endl;
//             }
//         }
//     }
    return triplets;
}

ConservativeSphericalPolygonInterpolation::Triplets ConservativeSphericalPolygonInterpolation::compute_2nd_order_triplets() {
    ATLAS_TRACE("ConservativeMethod::setup: build cons-2 interpolant matrix");
    ATLAS_ASSERT(not matrix_free_);
    const auto& src_points_ = data_->src_points_;
    const auto& src_iparam_ = data_->src_iparam_;
    Triplets triplets;
    size_t triplets_size    = 0;
    const auto& tgt_areas_v = data_->tgt_areas_;
    if (src_cell_data_) {
        const auto src_halo = array::make_view<int, 1>(src_mesh_.cells().halo());
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            triplets_size += (2 * nb_cells.size() + 1) * src_iparam_[scell].cell_idx.size();
        }
        triplets.reserve(triplets_size);
        std::vector<PointXYZ> Rsj;
        std::vector<PointXYZ> Aik;
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam  = src_iparam_[scell];
            if (iparam.cell_idx.size() == 0 && not src_halo(scell)) {
                continue;
            }
            /* // better conservation after Kritsikis et al. (2017)
            // NOTE: ommited here at cost of conservation due to more involved implementation in parallel
            // TEMP: use barycentre computed based on the source polygon's vertices, rather then
            //       the numerical barycentre based on barycentres of the intersections with target cells.
            // TODO: for a given source cell, collect the centroids of all its intersections with target cells
            //       to compute the numerical barycentre of the cell bases on intersection.
            PointXYZ Cs = {0., 0., 0.};
            for ( idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell ) {
                Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
            }
            */
            const PointXYZ& Cs = src_points_[scell];
            // compute gradient from cell
#if defined(__APPLE__)
volatile // On Apple, prevent FE_DIVBYZERO possibly triggered when inverting dual_area_inv a few lines below
#endif
            double dual_area_inv = 0.;
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            Rsj.resize(nb_cells.size());
            for (idx_t j = 0; j < nb_cells.size(); ++j) {
                idx_t nj         = next_index(j, nb_cells.size());
                idx_t sj         = nb_cells[j];
                idx_t nsj        = nb_cells[nj];
                const auto& Csj  = src_points_[sj];
                const auto& Cnsj = src_points_[nsj];
                if (ConvexSphericalPolygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16)) {
                    Rsj[j] = PointXYZ::cross(Cnsj, Csj);
                    dual_area_inv += ConvexSphericalPolygon({Cs, Csj, Cnsj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Csj, Cnsj);
                    dual_area_inv += ConvexSphericalPolygon({Cs, Cnsj, Csj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.; // dual_area_inv, even if protected, may still leads to FE_DIVBYZERO with Apple without volatile trick
            PointXYZ Rs   = {0., 0., 0.};
            for (idx_t j = 0; j < nb_cells.size(); ++j) {
                Rs = Rs + Rsj[j];
            }
            // assemble the matrix
            Aik.resize(iparam.cell_idx.size());
            for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                const PointXYZ& Csk   = iparam.centroids[icell];
                const PointXYZ Csk_Cs = Csk - Cs;
                Aik[icell]            = Csk_Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk_Cs));
                Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.tgt_weights[icell] * dual_area_inv);
            }
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    const idx_t tcell = iparam.cell_idx[icell];
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = next_index(j, nb_cells.size());
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        triplets.emplace_back(tcell, nsj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                    }
                    triplets.emplace_back(tcell, scell, iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                }
            }
            else {
                auto& tgt_csp2node_ = data_->tgt_csp2node_;
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    idx_t tcell            = iparam.cell_idx[icell];
                    idx_t tnode            = tgt_csp2node_[tcell];
                    double inv_node_weight = (tgt_areas_v[tnode] > 0.) ? 1. / tgt_areas_v[tnode] : 0.;
                    double csp2node_coef   = iparam.weights[icell] / iparam.tgt_weights[icell] * inv_node_weight;
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = next_index(j, nb_cells.size());
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        triplets.emplace_back(tnode, nsj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                    }
                    triplets.emplace_back(tnode, scell,
                                          (iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
                }
            }
        }
    }
    else {  // if ( not src_cell_data_ )
        auto& src_node2csp_ = data_->src_node2csp_;
        Workspace w;
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            const auto nb_nodes = get_node_neighbours(src_mesh_, snode, w);
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += (2 * nb_nodes.size() + 1) * src_iparam_[subcell].cell_idx.size();
            }
        }
        triplets.reserve(triplets_size);
        std::vector<PointXYZ> Rsj;
        std::vector<PointXYZ> Aik;
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            const auto nb_nodes = get_node_neighbours(src_mesh_, snode, w);
            // get the barycentre of the dual cell
            /* // better conservation
            PointXYZ Cs = {0., 0., 0.};
            for ( idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell ) {
                idx_t subcell      = src_node2csp_[snode][isubcell];
                const auto& iparam = src_iparam_[subcell];
                for ( idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell ) {
                    Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
                }
                const double Cs_norm = PointXYZ::norm( Cs );
                ATLAS_ASSERT( Cs_norm > 0. );
                Cs = PointXYZ::div( Cs, Cs_norm );
            }
            */
            const PointXYZ& Cs = src_points_[snode];
            // compute gradient from nodes
            double dual_area_inv = 0.;
            Rsj.resize(nb_nodes.size());
            const auto& Ns = src_points_[snode];
            for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                idx_t nj         = next_index(j, nb_nodes.size());
                idx_t sj         = nb_nodes[j];
                idx_t snj        = nb_nodes[nj];
                const auto& Nsj  = src_points_[sj];
                const auto& Nsnj = src_points_[snj];
                if (ConvexSphericalPolygon::GreatCircleSegment(Ns, Nsj).inLeftHemisphere(Nsnj, -1e-16)) {
                    Rsj[j] = PointXYZ::cross(Nsnj, Nsj);
                    dual_area_inv += ConvexSphericalPolygon({Ns, Nsj, Nsnj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Nsj, Nsnj);
                    dual_area_inv += ConvexSphericalPolygon({Ns, Nsnj, Nsj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.;
            PointXYZ Rs   = {0., 0., 0.};
            for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                Rs = Rs + Rsj[j];
            }
            // assemble the matrix
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell      = src_node2csp_[snode][isubcell];
                const auto& iparam = src_iparam_[subcell];
                if (iparam.cell_idx.size() == 0) {
                    continue;
                }
                Aik.resize(iparam.cell_idx.size());
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    const PointXYZ& Csk   = iparam.centroids[icell];
                    const PointXYZ Csk_Cs = Csk - Cs;
                    Aik[icell]            = Csk_Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk_Cs));
                    Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.tgt_weights[icell] * dual_area_inv);
                }
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                        const idx_t tcell = iparam.cell_idx[icell];
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = next_index(j, nb_nodes.size());
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                            triplets.emplace_back(tcell, snj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        }
                        triplets.emplace_back(tcell, snode, iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                    }
                }
                else {
                    auto& tgt_csp2node_ = data_->tgt_csp2node_;
                    for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                        idx_t tcell            = iparam.cell_idx[icell];
                        idx_t tnode            = tgt_csp2node_[tcell];
                        double inv_node_weight = (tgt_areas_v[tnode] > 1e-15) ? 1. / tgt_areas_v[tnode] : 0.;
                        double csp2node_coef = iparam.weights[icell] / iparam.tgt_weights[icell] * inv_node_weight;
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = next_index(j, nb_nodes.size());
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                            triplets.emplace_back(tnode, snj,
                                                  (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        }
                        triplets.emplace_back(
                            tnode, snode, (iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
                    }
                }
            }
        }
    }
    sort_and_accumulate_triplets(triplets);  // Very expensive!!! (90% of this routine). We need to avoid it
    return triplets;
}

void ConservativeSphericalPolygonInterpolation::do_execute(const FieldSet& src_fields, FieldSet& tgt_fields,
                                                           Metadata& metadata) const {
    std::vector<Metadata> md_array;
    md_array.resize( src_fields.size() );
    for (int i = 0; i < src_fields.size(); i++) { // TODO: memory-wise should we openmp this?
        do_execute(src_fields[i], tgt_fields[i], md_array[i]);
    }
    metadata = md_array[0]; // TODO: reduce metadata of a set of variables to a single metadata of a variable?
}

void ConservativeSphericalPolygonInterpolation::do_execute(const Field& src_field, Field& tgt_field,
                                                           Metadata& metadata) const {
    ATLAS_TRACE("ConservativeMethod::do_execute()");
    {
        if (src_field.dirty()) {
            ATLAS_TRACE("halo exchange source");
            src_field.haloExchange();
        }
    }
    StopWatch stopwatch;
    stopwatch.start();
    if (order_ == 1) {
        if (matrix_free_) {
            ATLAS_TRACE("matrix_free_order_1");
            const auto& src_iparam_ = data_->src_iparam_;
            const auto& tgt_areas_v = data_->tgt_areas_;

            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }
            const auto src_vals = array::make_view<double, 1>(src_field);
            auto tgt_vals       = array::make_view<double, 1>(tgt_field);
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                const auto& iparam = src_iparam_[scell];
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    tgt_vals(iparam.cell_idx[icell]) += iparam.weights[icell] * src_vals(scell);
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals[tcell] /= tgt_areas_v[tcell];
            }
        }
        else {
            ATLAS_TRACE("matrix_order_1");
            Method::do_execute(src_field, tgt_field, metadata);
        }
    }
    else if (order_ == 2) {
        if (matrix_free_) {
            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }
            ATLAS_TRACE("matrix_free_order_2");
            const auto& src_iparam_ = data_->src_iparam_;
            const auto& tgt_areas_v = data_->tgt_areas_;
            auto& src_points_ = data_->src_points_;
            const auto src_vals = array::make_view<double, 1>(src_field);
            auto tgt_vals       = array::make_view<double, 1>(tgt_field);
            const auto halo     = array::make_view<int, 1>(src_mesh_.cells().halo());
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                const auto& iparam       = src_iparam_[scell];
                if (iparam.cell_idx.size() == 0 or halo(scell)) {
                    continue;
                }
                const PointXYZ& Cs       = src_points_[scell];
                PointXYZ grad            = {0., 0., 0.};
                PointXYZ src_barycenter  = {0., 0., 0.};
                auto nb_cells = get_cell_neighbours(src_mesh_, scell);
                double dual_area_inv     = 0.;
                for (idx_t j = 0; j < nb_cells.size(); ++j) {
                    idx_t nj    = next_index(j, nb_cells.size());
                    idx_t sj     = nb_cells[j];
                    idx_t nsj    = nb_cells[nj];
                    const auto& Csj  = src_points_[sj];
                    const auto& Cnsj = src_points_[nsj];
                    double val = 0.5 * (src_vals(nj) + src_vals(nsj)) - src_vals(j);
                    if (ConvexSphericalPolygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16)) {
                        dual_area_inv += ConvexSphericalPolygon({Cs, Csj, Cnsj}).area();
                    }
                    else {
                        //val *=  -1;
                        dual_area_inv += ConvexSphericalPolygon({Cs, Cnsj, Csj}).area();
                    }
                    grad = grad + PointXYZ::mul(PointXYZ::cross(Csj, Cnsj), val);
                }
                dual_area_inv = ((dual_area_inv > 0.) ? 1. / dual_area_inv : 0.);
                grad = PointXYZ::mul(grad, dual_area_inv);
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    src_barycenter = src_barycenter + PointXYZ::mul(iparam.centroids[icell], iparam.weights[icell]);
                }
                src_barycenter = PointXYZ::div(src_barycenter, PointXYZ::norm(src_barycenter));
                grad           = grad - PointXYZ::mul(src_barycenter, PointXYZ::dot(grad, src_barycenter));
                ATLAS_ASSERT(std::abs(PointXYZ::dot(grad, src_barycenter)) < 1e-14);
                for (idx_t icell = 0; icell < iparam.cell_idx.size(); ++icell) {
                    tgt_vals(iparam.cell_idx[icell]) +=
                        iparam.weights[icell] *
                        (src_vals(scell) + PointXYZ::dot(grad, iparam.centroids[icell] - src_barycenter));
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals[tcell] /= tgt_areas_v[tcell];
            }
        }
        else {
            ATLAS_TRACE("matrix_order_2");
            Method::do_execute(src_field, tgt_field, metadata);
        }
    }

    stopwatch.stop();
    {
        ATLAS_TRACE("halo exchange target");
        if (tgt_field.hostNeedsUpdate()) {
            tgt_field.updateHost();
        }
        tgt_field.haloExchange();
        tgt_field.setDeviceNeedsUpdate(true);
    }

    auto remap_stat = remap_stat_;
    if (statistics_conservation_) {
        const auto src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
        const auto src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
        const auto src_node_halo  = array::make_view<int, 1>(src_mesh_.nodes().halo());
        const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
        const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
        const auto tgt_node_halo  = array::make_view<int, 1>(tgt_mesh_.nodes().halo());
        const auto& src_areas_v   = data_->src_areas_;
        const auto& tgt_areas_v   = data_->tgt_areas_;

        const auto src_vals = array::make_view<double, 1>(src_field);
        const auto tgt_vals = array::make_view<double, 1>(tgt_field);

        double err_remap_cons     = 0.;
        if (src_cell_data_) {
            for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
                if (src_cell_halo(spt)) {
                    continue;
                }
                err_remap_cons += src_vals(spt) * src_areas_v[spt];
            }
        }
        else {
            for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
                if (src_node_halo(spt) or src_node_ghost(spt)) {
                    continue;
                }
                err_remap_cons += src_vals(spt) * src_areas_v[spt];
            }
        }
        if (tgt_cell_data_) {
            for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
                if (tgt_cell_halo(tpt)) {
                    continue;
                }
                err_remap_cons -= tgt_vals(tpt) * tgt_areas_v[tpt];
            }
        }
        else {
            for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
                if (tgt_node_halo(tpt) or tgt_node_ghost(tpt)) {
                    continue;
                }
                err_remap_cons -= tgt_vals(tpt) * tgt_areas_v[tpt];
            }
        }
        ATLAS_TRACE_MPI(ALLREDUCE) { mpi::comm().allReduceInPlace(&err_remap_cons, 1, eckit::mpi::sum()); }
        remap_stat.errors[Statistics::Errors::REMAP_CONS] = err_remap_cons / unit_sphere_area();

        metadata.set("conservation_error", remap_stat.errors[Statistics::Errors::REMAP_CONS]);
    }
    if (statistics_intersection_) {
        metadata.set("polygons.source", remap_stat.counts[Statistics::Counts::SRC_PLG]);
        metadata.set("polygons.target", remap_stat.counts[Statistics::Counts::TGT_PLG]);
        metadata.set("polygons.intersections", remap_stat.counts[Statistics::Counts::INT_PLG]);
        metadata.set("polygons.uncovered_source", remap_stat.counts[Statistics::UNCVR_SRC]);
        if (validate_) {
            metadata.set("source_area_error.L1", remap_stat.errors[Statistics::Errors::SRC_INTERSECTPLG_L1]);
            metadata.set("source_area_error.Linf", remap_stat.errors[Statistics::Errors::SRC_INTERSECTPLG_LINF]);
            metadata.set("target_area_error.L1", remap_stat.errors[Statistics::Errors::TGT_INTERSECTPLG_L1]);
            metadata.set("target_area_error.Linf", remap_stat.errors[Statistics::Errors::TGT_INTERSECTPLG_LINF]);
        }
    }

    if (statistics_intersection_ || statistics_conservation_) {
        remap_stat.fillMetadata(metadata);
    }

    auto& timings = data_->timings;
    if (statistics_timings_) {
        metadata.set("timings.target_kdtree_search", timings.target_kdtree_search);
        metadata.set("timings.polygon_intersections", timings.polygon_intersections);
    }
    metadata.set("timings.source_polygons_assembly", timings.source_polygons_assembly);
    metadata.set("timings.target_polygons_assembly", timings.target_polygons_assembly);
    metadata.set("timings.target_kdtree_assembly", timings.target_kdtree_assembly);
    metadata.set("timings.source_polygons_filter", timings.source_polygons_filter);
    metadata.set("timings.matrix_assembly", timings.matrix_assembly);
    metadata.set("timings.interpolation", stopwatch.elapsed());

    metadata.set("memory.matrix", matrix_free_ ? 0 : matrix().footprint());
    metadata.set("memory.src_points", memory_of(data_->src_points_));
    metadata.set("memory.tgt_points", memory_of(data_->tgt_points_));
    metadata.set("memory.src_areas", memory_of(data_->src_points_));
    metadata.set("memory.tgt_areas", memory_of(data_->tgt_areas_));
    metadata.set("memory.src_csp2node", memory_of(data_->src_csp2node_));
    metadata.set("memory.tgt_csp2node", memory_of(data_->tgt_csp2node_));
    metadata.set("memory.src_node2csp", memory_of(data_->src_node2csp_));
    metadata.set("memory.tgt_node2csp", memory_of(data_->tgt_node2csp_));
    metadata.set("memory.src_iparam", memory_of(data_->src_iparam_));
}

void ConservativeSphericalPolygonInterpolation::print(std::ostream& out) const {
    out << "ConservativeMethod{";
    out << "order:" << order_;
    int halo = 0;
    src_mesh_.metadata().get("halo", halo);
    out << ", source:" << (src_cell_data_ ? "cells(" : "nodes(") << src_mesh_.grid().name() << ",halo=" << halo << ")";
    tgt_mesh_.metadata().get("halo", halo);
    out << ", target:" << (tgt_cell_data_ ? "cells(" : "nodes(") << tgt_mesh_.grid().name() << ",halo=" << halo << ")";
    out << ", normalise_intersections:" << normalise_intersections_;
    out << ", matrix_free:" << matrix_free_;
    out << ", statistics.intersection:" << statistics_intersection_;
    out << ", statistics.conservation:" << statistics_conservation_;
    out << ", cached_matrix:" << not(matrixAllocated() || matrix_free_);
    out << ", cached_data:" << bool(sharable_data_.use_count() == 0);
    size_t footprint{};
    if (not matrix_free_) {
        footprint += matrix().footprint();
    }
    footprint += data_->footprint();
    out << ", footprint:" << eckit::Bytes(footprint);
    out << "}";
}

Cache ConservativeSphericalPolygonInterpolation::createCache() const {
    interpolation::Cache cache;
    if (not matrix_free_) {
        cache.add(Method::createCache());
    }
    cache.add(cache_);
    return cache;
}

void ConservativeSphericalPolygonInterpolation::setup_stat() const {
    const auto src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto& src_areas_v   = data_->src_areas_;
    const auto& tgt_areas_v   = data_->tgt_areas_;
    double geo_create_err     = 0.;
    double src_tgt_sums[2]    = {0., 0.};
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_areas_v.size(); ++spt) {
            if (not src_cell_halo(spt)) {
                src_tgt_sums[0] += src_areas_v[spt];
            }
        }
    }
    else {
        for (idx_t src = 0; src < src_areas_v.size(); ++src) {
            if (not src_node_ghost(src)) {
                src_tgt_sums[0] += src_areas_v[src];
            }
        }
    }
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_cell_halo(tpt)) {
                src_tgt_sums[1] += tgt_areas_v[tpt];
            }
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_node_ghost(tpt)) {
                src_tgt_sums[1] += tgt_areas_v[tpt];
            }
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) { mpi::comm().allReduceInPlace(src_tgt_sums, 2, eckit::mpi::sum()); }

    remap_stat_.src_area_sum = src_tgt_sums[0];
    remap_stat_.tgt_area_sum = src_tgt_sums[1];

    geo_create_err                                   = std::abs(src_tgt_sums[0] - src_tgt_sums[1]) / unit_sphere_area();
    remap_stat_.errors[Statistics::Errors::SRCTGT_INTERSECTPLG_DIFF] = geo_create_err;
}

Field ConservativeSphericalPolygonInterpolation::Statistics::diff(const Interpolation& interpolation,
                                                                  const Field source, const Field target) {
    Field diff     = interpolation.source().createField(source, option::name("diff"));
    auto diff_vals = array::make_view<double, 1>(diff);

    const auto src_vals = array::make_view<double, 1>(source);
    const auto tgt_vals = array::make_view<double, 1>(target);

    auto cachable_data_       = ConservativeSphericalPolygonInterpolation::Cache(interpolation).get();
    const auto& src_areas_v   = cachable_data_->src_areas_;
    const auto& tgt_csp2node_ = cachable_data_->tgt_csp2node_;
    const auto& src_node2csp_ = cachable_data_->src_node2csp_;
    const auto& src_iparam_   = cachable_data_->src_iparam_;
    const auto& src_mesh_     = extract_mesh(cachable_data_->src_fs_);
    const auto& tgt_mesh_     = extract_mesh(cachable_data_->tgt_fs_);
    const auto src_cell_data_ = bool(functionspace::CellColumns(interpolation.source()));
    const auto tgt_cell_data_ = bool(functionspace::CellColumns(interpolation.target()));
    const auto src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto src_node_halo  = array::make_view<int, 1>(src_mesh_.nodes().halo());
    const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_cell_halo(spt)) {
                continue;
            }
            double diff        = src_vals(spt) * src_areas_v[spt];
            const auto& iparam = src_iparam_[spt];
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    if (tgt_cell_halo(tcell) < 1) {
                        diff -= tgt_vals(tcell) * iparam.weights[icell];
                    }
                }
            }
            else {
                for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    idx_t tnode = tgt_csp2node_[tcell];
                    if (not tgt_node_ghost(tnode)) {
                        diff -= tgt_vals(tnode) * iparam.weights[icell];
                    }
                }
            }
            if (src_areas_v[spt] > 0.) {
                diff_vals(spt) = diff / src_areas_v[spt];
            }
            else {
                Log::info() << " at cell " << spt << " cell-area: " << src_areas_v[spt] << std::endl;
                diff_vals(spt) = std::numeric_limits<double>::max();
            }
        }
    }
    else {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_node_ghost(spt) or src_node_halo(spt) != 0) {
                continue;
            }
            double diff          = src_vals(spt) * src_areas_v[spt];
            const auto& node2csp = src_node2csp_[spt];
            for (idx_t subcell = 0; subcell < node2csp.size(); ++subcell) {
                const auto& iparam = src_iparam_[node2csp[subcell]];
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                        idx_t tcell = iparam.cell_idx[icell];
                        if (tgt_cell_halo(tcell) < 1) {
                            diff -= tgt_vals(tcell) * iparam.weights[icell];
                        }
                    }
                }
                else {
                    for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                        idx_t tcell = iparam.cell_idx[icell];
                        idx_t tnode = tgt_csp2node_[tcell];
                        if (not tgt_node_ghost(tnode)) {
                            diff -= tgt_vals(tnode) * iparam.weights[icell];
                        }
                    }
                }
            }
            if (src_areas_v[spt] > 0.) {
                diff_vals(spt) = diff / src_areas_v[spt];
            }
            else {
                diff_vals(spt) = std::numeric_limits<double>::max();
                Log::info() << " at cell " << spt << " cell-area: " << src_areas_v[spt] << std::endl;
            }
        }
    }
    return diff;
}


ConservativeSphericalPolygonInterpolation::Metadata ConservativeSphericalPolygonInterpolation::Statistics::accuracy(const Interpolation& interpolation,
                                                                         const Field target,
                                                                         std::function<double(const PointLonLat&)> func) {
    auto tgt_vals             = array::make_view<double, 1>(target);
    auto cachable_data_       = ConservativeSphericalPolygonInterpolation::Cache(interpolation).get();
    auto tgt_mesh_            = extract_mesh(cachable_data_->src_fs_);
    auto tgt_cell_data_       = extract_mesh(cachable_data_->tgt_fs_);
    const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    const auto tgt_node_halo  = array::make_view<int, 1>(tgt_mesh_.nodes().halo());
    const auto& tgt_areas_v   = cachable_data_->tgt_areas_;
    auto& tgt_points_         = cachable_data_->tgt_points_;
    double err_remap_l2       = 0.;
    double err_remap_linf     = 0.;
    if (tgt_cell_data_) {
        size_t ncells = std::min<size_t>(tgt_vals.size(), tgt_mesh_.cells().size());
        for (idx_t tpt = 0; tpt < ncells; ++tpt) {
            ATLAS_ASSERT(tpt < tgt_cell_halo.size());
            if (tgt_cell_halo(tpt)) {
                continue;
            }
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas_v[tpt];
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    else {
        size_t nnodes = std::min<size_t>(tgt_vals.size(), tgt_mesh_.nodes().size());
        for (idx_t tpt = 0; tpt < nnodes; ++tpt) {
            if (tgt_node_ghost(tpt) or tgt_node_halo(tpt)) {
                continue;
            }
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas_v[tpt];
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&err_remap_l2, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&err_remap_linf, 1, eckit::mpi::max());
    }
    errors[Statistics::Errors::REMAP_L2]   = std::sqrt(err_remap_l2 / unit_sphere_area());
    errors[Statistics::Errors::REMAP_LINF] = err_remap_linf;
    Metadata metadata;
    metadata.set("errors.REMAP_L2", errors[Statistics::Errors::REMAP_L2]);
    metadata.set("errors.REMAP_LINF", errors[Statistics::Errors::REMAP_LINF]);
    return metadata;
}

void ConservativeSphericalPolygonInterpolation::dump_intersection(const std::string msg,
                                                                  const ConvexSphericalPolygon& plg_1,
                                                                  const CSPolygonArray& plg_2_array,
                                                                  const std::vector<idx_t>& plg_2_idx_array) const {
#if PRINT_BAD_POLYGONS
    double plg_1_coverage = 0.;
    std::vector<double> int_areas;
    int_areas.resize(plg_2_idx_array.size());
    for (int i = 0; i < plg_2_idx_array.size(); ++i) {
        const auto plg_2_idx = plg_2_idx_array[i];
        const auto& plg_2    = std::get<0>(plg_2_array[plg_2_idx]);
        auto iplg            = plg_1.intersect(plg_2, false, 1.e5 * std::numeric_limits<double>::epsilon());
        plg_1_coverage += iplg.area();
        int_areas[i] = iplg.area();
    }
    if (std::abs(plg_1.area() - plg_1_coverage) > 0.01 * plg_1.area()) {
        Log::info().indent();
        Log::info() << msg << ", Polygon_1.area: " << plg_1.area() << ", covered: " << plg_1_coverage << std::endl;
        Log::info() << "Polygon_1 : ";
        Log::info().precision(18);
        plg_1.print(Log::info());
        Log::info() << std::endl << "Printing " << plg_2_idx_array.size() << " covering polygons -->" << std::endl;
        Log::info().indent();
        plg_1_coverage = plg_1.area();
        for (int i = 0; i < plg_2_idx_array.size(); ++i) {
            const auto plg_2_idx = plg_2_idx_array[i];
            const auto& plg_2    = std::get<0>(plg_2_array[plg_2_idx]);
            plg_2.print(Log::info());
            Log::info() << "\ninters:";
            plg_1_coverage -= int_areas[i];
            auto iplg       = plg_1.intersect(plg_2, false, 1.e5 * std::numeric_limits<double>::epsilon());
            Log::info().indent();
            iplg.print(Log::info());
            Log::info() << ", inters.-area: " << iplg.area() << ", plg1-area left: " << plg_1_coverage << std::endl;
            Log::info().unindent();
            Log::info() << std::endl;
        }
        Log::info() << std::endl;
        Log::info().unindent();
        Log::info().unindent();
    }
#endif
}

template <class TargetCellsIDs>
void ConservativeSphericalPolygonInterpolation::dump_intersection(const std::string msg,
                                                                  const ConvexSphericalPolygon& plg_1,
                                                                  const CSPolygonArray& plg_2_array,
                                                                  const TargetCellsIDs& plg_2_idx_array) const {
#if PRINT_BAD_POLYGONS
    std::vector<idx_t> idx_array;
    idx_array.resize(plg_2_idx_array.size());
    for (int i = 0; i < plg_2_idx_array.size(); ++i) {
        idx_array[i] = plg_2_idx_array[i].payload();
    }
    dump_intersection(msg, plg_1, plg_2_array, idx_array);
#endif
}

ConservativeSphericalPolygonInterpolation::Cache::Cache(std::shared_ptr<InterpolationCacheEntry> entry):
    interpolation::Cache(entry), entry_(dynamic_cast<Data*>(entry.get())) {}

ConservativeSphericalPolygonInterpolation::Cache::Cache(const interpolation::Cache& c):
    interpolation::Cache(c, Data::static_type()), entry_{dynamic_cast<const Data*>(c.get(Data::static_type()))} {}

ConservativeSphericalPolygonInterpolation::Cache::Cache(const Interpolation& interpolation):
    Cache(interpolation::Cache(interpolation)) {}

size_t ConservativeSphericalPolygonInterpolation::Data::footprint() const {
    size_t mem_total{0};
    mem_total += memory_of(src_points_);
    mem_total += memory_of(tgt_points_);
    mem_total += memory_of(src_areas_);
    mem_total += memory_of(tgt_areas_);
    mem_total += memory_of(src_csp2node_);
    mem_total += memory_of(tgt_csp2node_);
    mem_total += memory_of(src_node2csp_);
    mem_total += memory_of(tgt_node2csp_);
    mem_total += memory_of(src_iparam_);
    return mem_total;
}


void ConservativeSphericalPolygonInterpolation::Data::print(std::ostream& out) const {
    out << "Memory usage of ConservativeMethod: " << eckit::Bytes(footprint()) << "\n";
    out << "- src_points_   \t" << eckit::Bytes(memory_of(src_points_)) << "\n";
    out << "- tgt_points_   \t" << eckit::Bytes(memory_of(tgt_points_)) << "\n";
    out << "- src_areas_    \t" << eckit::Bytes(memory_of(src_areas_)) << "\n";
    out << "- tgt_areas_    \t" << eckit::Bytes(memory_of(tgt_areas_)) << "\n";
    out << "- src_csp2node_ \t" << eckit::Bytes(memory_of(src_csp2node_)) << "\n";
    out << "- tgt_csp2node_ \t" << eckit::Bytes(memory_of(tgt_csp2node_)) << "\n";
    out << "- src_node2csp_ \t" << eckit::Bytes(memory_of(src_node2csp_)) << "\n";
    out << "- tgt_node2csp_ \t" << eckit::Bytes(memory_of(tgt_node2csp_)) << "\n";
    out << "- src_iparam_   \t" << eckit::Bytes(memory_of(src_iparam_)) << "\n";
}

void ConservativeSphericalPolygonInterpolation::Statistics::fillMetadata(Metadata& metadata) {
    // errors
    metadata.set("errors.SRC_SUBPLG_L1", errors[SRC_SUBPLG_L1]);
    metadata.set("errors.SRC_SUBPLG_LINF", errors[SRC_SUBPLG_LINF]);
    metadata.set("errors.TGT_SUBPLG_L1", errors[TGT_SUBPLG_L1]);
    metadata.set("errors.TGT_SUBPLG_LINF", errors[TGT_SUBPLG_LINF]);
    metadata.set("errors.SRC_INTERSECTPLG_L1", errors[SRC_INTERSECTPLG_L1]);
    metadata.set("errors.SRC_INTERSECTPLG_LINF", errors[SRC_INTERSECTPLG_LINF]);
    metadata.set("errors.TGT_INTERSECTPLG_L1", errors[TGT_INTERSECTPLG_L1]);
    metadata.set("errors.TGT_INTERSECTPLG_LINF", errors[TGT_INTERSECTPLG_LINF]);
    metadata.set("errors.SRCTGT_INTERSECTPLG_DIFF", errors[SRCTGT_INTERSECTPLG_DIFF]);
    metadata.set("errors.REMAP_CONS", errors[REMAP_CONS]);
    metadata.set("errors.REMAP_L2", errors[REMAP_L2]);
    metadata.set("errors.REMAP_LINF", errors[REMAP_LINF]);
}

ConservativeSphericalPolygonInterpolation::Statistics::Statistics() {
    std::fill(std::begin(counts), std::end(counts), 0);
    std::fill(std::begin(errors), std::end(errors), 0.);
}

ConservativeSphericalPolygonInterpolation::Statistics::Statistics(const Metadata& metadata): Statistics() {
    // errors
    metadata.get("errors.SRC_SUBPLG_L1", errors[SRC_SUBPLG_L1]);
    metadata.get("errors.SRC_SUBPLG_LINF", errors[SRC_SUBPLG_LINF]);
    metadata.get("errors.TGT_SUBPLG_L1", errors[TGT_SUBPLG_L1]);
    metadata.get("errors.TGT_SUBPLG_LINF", errors[TGT_SUBPLG_LINF]);
    metadata.get("errors.SRC_INTERSECTPLG_L1", errors[SRC_INTERSECTPLG_L1]);
    metadata.get("errors.SRC_INTERSECTPLG_LINF", errors[SRC_INTERSECTPLG_LINF]);
    metadata.get("errors.TGT_INTERSECTPLG_L1", errors[TGT_INTERSECTPLG_L1]);
    metadata.get("errors.TGT_INTERSECTPLG_LINF", errors[TGT_INTERSECTPLG_LINF]);
    metadata.get("errors.SRCTGT_INTERSECTPLG_DIFF", errors[SRCTGT_INTERSECTPLG_DIFF]);
    metadata.get("errors.REMAP_CONS", errors[REMAP_CONS]);
    metadata.get("errors.REMAP_L2", errors[REMAP_L2]);
    metadata.get("errors.REMAP_LINF", errors[REMAP_LINF]);
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
