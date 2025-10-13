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
#include <unordered_set>
#include <vector>

#include "ConservativeSphericalPolygonInterpolation.h"

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
#include "atlas/util/detail/filesystem.h"

#include "eckit/log/Bytes.h"
#include "eckit/log/ProgressTimer.h"

#define PRINT_BAD_POLYGONS 0

namespace atlas {
namespace interpolation {
namespace method {

MethodBuilder<ConservativeSphericalPolygonInterpolation> __builder("conservative-spherical-polygon");

using runtime::trace::StopWatch;
using Polygon = util::ConvexSphericalPolygon;
using PolygonArray = std::vector<util::ConvexSphericalPolygon>;

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


template<typename PolygonContainer>
std::string to_json(const PolygonContainer& polygons, int precision = 0) {
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


template<typename Polygons>
std::string polygons_to_json(const Polygons& polygons, int precision = 0) {
    return polygons_to_json(polygons.begin(),polygons.end(),precision);
}


template<typename Polygons>
void dump_polygons_to_json( const Polygon& t_csp,
                    double pointsSameEPS,
                    const Polygons& source_polygons,
                    const std::vector<idx_t>& source_polygons_considered_indices,
                    const std::string folder,
                    const std::string name) {
    PolygonArray csp_arr{ t_csp };
    PolygonArray csp_arr_intersecting {t_csp};
    PolygonArray intersections;
    int count = 1;
    for( auto& s_idx : source_polygons_considered_indices ) {
        auto s_csp = source_polygons[s_idx];
        csp_arr.emplace_back( s_csp );
        std::fstream file_plg_debug(folder + name + "_" + std::to_string(count++) + ".debug", std::ios::out);
        Polygon iplg = t_csp.intersect(s_csp, &file_plg_debug, pointsSameEPS);
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
        mem += memory_of(params.csp_ids);
        mem += memory_of(params.centroids);
        mem += memory_of(params.weights);
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


#if ATLAS_BUILD_TYPE_DEBUG
int inside_vertices(const Polygon& plg1, const Polygon& plg2, int& pout) {
    int points_in = 0;
    pout          = 0;
    for (int j = 0; j < plg2.size(); j++) {
        int i = 0;
        for (; i < plg1.size(); i++) {
            int in   = (i != plg1.size() - 1) ? i + 1 : 0;
            auto gss = Polygon::GreatCircleSegment(plg1[i], plg1[in]);
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
#endif


inline bool valid_point(idx_t node_idx, const array::ArrayView<int, 1>& node_flags) {
    return not util::Bitflags::view(node_flags(node_idx)).check(util::Topology::INVALID);
}

}  // namespace


ConservativeSphericalPolygonInterpolation::ConservativeSphericalPolygonInterpolation(const Config& config):
    Method(config), validate_(false), src_cell_data_(true), tgt_cell_data_(true), normalise_(false), limit_("none"),
    order_(1), matrix_free_(false), n_spoints_(0), n_tpoints_(0) {
    config.get("validate", validate_ = false);
    config.get("order", order_ = 1);
    config.get("normalise", normalise_ = false);
    config.get("limit", limit_ = "none");
    config.get("matrix_free", matrix_free_ = false);
    config.get("src_cell_data", src_cell_data_ = true);
    config.get("tgt_cell_data", tgt_cell_data_ = true);
    config.get("statistics.all", remap_stat_.all = false);
    config.get("statistics.accuracy", remap_stat_.accuracy = false);
    config.get("statistics.conservation", remap_stat_.conservation = false);
    config.get("statistics.intersection", remap_stat_.intersection = false);
    config.get("statistics.timings", remap_stat_.timings = false);
    if (remap_stat_.all) {
        Log::warning() << "statistics.all required. Enabling validate, statistics.timings, statistics.intersection, statistics.conservation, and statistics.accuracy." << std::endl;
        validate_ = true;
        remap_stat_.accuracy = true;
        remap_stat_.conservation = true;
        remap_stat_.intersection = true;
        remap_stat_.timings = true;
    }
    if (remap_stat_.intersection) {
        Log::warning() << "statistics.intersection required. Enabling validate." << std::endl;
        validate_ = true;
    }

    sharable_data_ = std::make_shared<Data>();
    cache_         = Cache(sharable_data_);
    data_          = sharable_data_.get();

    ATLAS_ASSERT(sharable_data_.use_count() == 2);
}


inline int ConservativeSphericalPolygonInterpolation::next_index(int current_index, int size) const {
#if ATLAS_BUILD_TYPE_DEBUG
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
#endif
    return (current_index < size - 1) ? current_index + 1 : current_index + 1 - size;
}


inline int ConservativeSphericalPolygonInterpolation::prev_index(int current_index, int size) const {
#if ATLAS_BUILD_TYPE_DEBUG
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
#endif
    return (current_index >= 1) ? current_index - 1 : current_index - 1 + size;
}


struct ConservativeSphericalPolygonInterpolation::Workspace_get_cell_neighbours {
    PointXYZ p0;
    PointLonLat p0_ll;
    PointXYZ p1;
    PointLonLat p1_ll;
};


// get counter-clockwise sorted neighbours of a cell
std::vector<idx_t> ConservativeSphericalPolygonInterpolation::get_cell_neighbours(Mesh& mesh, idx_t cell, Workspace_get_cell_neighbours& w) const {
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

    auto& p0 = w.p0;
    auto& p1 = w.p1;
    auto& p0_ll = w.p0_ll;
    auto& p1_ll = w.p1_ll;

    for (idx_t inode = 0; inode < n_nodes; ++inode) {
        idx_t node0             = cell2node(cell, inode);
        idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
        p0_ll[0] = nodes_ll(node0, 0);
        p0_ll[1] = nodes_ll(node0, 1);
        p1_ll[0] = nodes_ll(node1, 0);
        p1_ll[1] = nodes_ll(node1, 1);
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


struct ConservativeSphericalPolygonInterpolation::Workspace_get_node_neighbours {
    std::vector<idx_t> nbr_nodes_od;
    std::vector< std::array<idx_t,2> > cnodes;
};


// get cyclically sorted node neighbours without using edge connectivity
std::vector<idx_t> ConservativeSphericalPolygonInterpolation::get_node_neighbours(Mesh& mesh, idx_t node_id, Workspace_get_node_neighbours& w) const {
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


void ConservativeSphericalPolygonInterpolation::
init_polygons_data(
    // input
    FunctionSpace fs,
    // output
    Data::PolygonsData& md) {
    using Context = Data::PolygonsData::Context;

    md.cell_data = functionspace::CellColumns(fs);
    std::vector<idx_t>& csp2node = md.csp2node;
    std::vector<std::vector<idx_t>>& node2csp = md.node2csp;
    gidx_t& csp_size = md.csp_size;
    std::vector<idx_t>& csp_cell_index = md.csp_cell_index;
    std::vector<idx_t>& csp_index = md.csp_index;

    auto mesh = extract_mesh(fs);
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto cell_halo   = array::make_view<int, 1>(mesh.cells().halo());
    const auto cell_flags  = array::make_view<int, 1>(mesh.cells().flags());
    const auto cell_invalid = [&cell_flags](idx_t cell) {
        return util::Bitflags::view(cell_flags(cell)).check(util::Topology::INVALID);
    };
    auto fs_halo = md.cell_data ? functionspace::CellColumns(fs).halo().size() : functionspace::NodeColumns(fs).halo().size();
    csp_size = 0;
    csp_cell_index.reserve(mesh.cells().size());
    csp_index.reserve(mesh.cells().size());
    if (md.cell_data) {
        auto skip_cell = [&](idx_t cell) {
            if (md.context == Context::SOURCE) {
                return cell_halo(cell) > fs_halo || cell_invalid(cell);
            }
            else {
                return cell_halo(cell) > 0;
            }
        };
        for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
            if (skip_cell(cell)) {
                continue;
            }
            csp_cell_index.emplace_back(cell);
            csp_index.emplace_back(cell);
            csp_size++;
        }
    }
    else {
        auto skip_cell = [&](idx_t cell) {
            if (md.context == Context::SOURCE) {
                return cell_halo(cell) > fs_halo || cell_invalid(cell);
            }
            else {
                return cell_halo(cell) > 1 || cell_invalid(cell);
            }
        };

        const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
        const auto node_flags = array::make_view<int, 1>(mesh.nodes().flags());
        auto ll2xyz = [](const atlas::PointLonLat& p_ll) {
            PointXYZ p_xyz;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p_ll, p_xyz);
            return p_xyz;
        };
        csp_index.emplace_back(0);
        for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
            if (skip_cell(cell)) {
                continue;
            }
            const idx_t n_nodes = cell2node.cols(cell);
            for (idx_t inode = 0; inode < n_nodes; ++inode) {
                idx_t node0             = cell2node(cell, inode);
                if (not valid_point(node0, node_flags)) {
                    continue;
                }
                idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
                const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
                const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
                PointXYZ p0             = ll2xyz(p0_ll);
                PointXYZ p1             = ll2xyz(p1_ll);
                if (PointXYZ::norm(p0 - p1) < 1e-14) {
                    continue;  // skip this edge = a pole point
                }
                csp_size++;
            }
            csp_index.emplace_back(csp_size);
            csp_cell_index.emplace_back(cell);
        }
        csp2node.resize(csp_size);
        std::fill(csp2node.begin(), csp2node.end(), -1);
        node2csp.resize(mesh.nodes().size());
    }
}


ConservativeSphericalPolygonInterpolation::Polygon ConservativeSphericalPolygonInterpolation::
get_csp_celldata(idx_t csp_id, const Mesh& mesh, const Data::PolygonsData& md) {

    // TODO: Performance optimisations:
    //    - Do not extract array views, mesh connectivity for each csp_id
    //    - pts_ll should also not be created/destroyed for each csp_id

    const auto& cell2node  = mesh.cells().node_connectivity();
    std::vector<PointLonLat> pts_ll;

    idx_t cell = csp_to_cell(csp_id, md);

    const idx_t n_nodes = cell2node.cols(cell);
    const auto nodes_ll   = array::make_view<const double, 2>(mesh.nodes().lonlat());
    pts_ll.clear();
    pts_ll.resize(n_nodes);

    for (idx_t jnode = 0; jnode < n_nodes; ++jnode) {
        idx_t inode   = cell2node(cell, jnode);
        pts_ll[jnode] = PointLonLat{nodes_ll(inode, 0), nodes_ll(inode, 1)};
    }
    return Polygon(pts_ll);
}


ConservativeSphericalPolygonInterpolation::Polygon ConservativeSphericalPolygonInterpolation::
get_csp_nodedata(idx_t csp_id, const Mesh& mesh, Data::PolygonsData& md ) {
    std::vector<idx_t>& csp2node = md.csp2node;
    std::vector<std::vector<idx_t>>& node2csp = md.node2csp;

    // TODO: Performance optimisations:
    //    - Do not extract array views, mesh connectivity for each csp_id
    //    - pts_ll, pts_xyz, pts_idx, subpol_pts_ll should also not be created/destroyed for each csp_id

    const auto& cell2node  = mesh.cells().node_connectivity();
    const auto node_flags = array::make_view<int, 1>(mesh.nodes().flags());
    std::vector<PointLonLat> pts_ll;

    idx_t cell, isubcell;
    std::tie(cell,isubcell) = csp_to_cell_and_subcell(csp_id, md);

    // idx_t cell;
    // auto isubcell = get_owner_cell_nodedata(csp_id, csp_cell_index, csp_index, cell);

    const idx_t n_nodes = cell2node.cols(cell);
    const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    pts_ll.clear();
    pts_ll.resize(n_nodes);

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

    PointXYZ cell_mid(0., 0., 0.);  // cell centre
    std::vector<PointXYZ> pts_xyz;
    std::vector<int> pts_idx;
#if ATLAS_BUILD_TYPE_DEBUG
    ATLAS_ASSERT(cell < cell2node.rows());
    ATLAS_ASSERT(n_nodes > 2);
#endif
    pts_xyz.clear();
    pts_ll.clear();
    pts_idx.clear();
    pts_xyz.reserve(n_nodes);
    pts_ll.reserve(n_nodes);
    pts_idx.reserve(n_nodes);
    for (idx_t inode = 0; inode < n_nodes; ++inode) {
        idx_t node0             = cell2node(cell, inode);
        if (not valid_point(node0, node_flags)) {
            continue;
        }
        idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
        const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
        const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
        PointXYZ p0             = ll2xyz(p0_ll);
        PointXYZ p1             = ll2xyz(p1_ll);
        if (PointXYZ::norm(p0 - p1) < 1e-14) {
            continue;  // skip this edge = a pole point
        }
        pts_xyz.emplace_back(p0);
        pts_ll.emplace_back(p0_ll);
        pts_idx.emplace_back(inode);
        cell_mid = cell_mid + p0;
        cell_mid = cell_mid + p1;
    }
    cell_mid  = PointXYZ::normalize(cell_mid);
#if ATLAS_BUILD_TYPE_DEBUG
    ATLAS_ASSERT(pts_xyz.size() > 2);
    ATLAS_ASSERT(isubcell < pts_idx.size());
#endif
    PointLonLat cell_ll = xyz2ll(cell_mid);
    idx_t inode = isubcell;
    int inode_n        = next_index(inode, pts_idx.size());
    PointXYZ iedge_mid = PointXYZ::normalize(pts_xyz[inode] + pts_xyz[inode_n]);
    int inode_nn = next_index(inode_n, pts_idx.size());
    if (PointXYZ::norm(pts_xyz[inode_nn] - pts_xyz[inode_n]) < 1e-14) {
        ATLAS_THROW_EXCEPTION("Three cell vertices on a same great arc!");
    }
    PointXYZ jedge_mid;
    jedge_mid = pts_xyz[inode_nn] + pts_xyz[inode_n];
    jedge_mid = PointXYZ::div(jedge_mid, PointXYZ::norm(jedge_mid));
    std::array<PointLonLat, 4> subpol_pts_ll;
    subpol_pts_ll[0] = cell_ll;
    subpol_pts_ll[1] = xyz2ll(iedge_mid);
    subpol_pts_ll[2] = pts_ll[inode_n];
    subpol_pts_ll[3] = xyz2ll(jedge_mid);

    idx_t node_n       = cell2node(cell, inode_n);
    if (csp2node[csp_id] == -1) {
        csp2node[csp_id] = node_n;
        node2csp[node_n].emplace_back(csp_id);
    }
    return Polygon(subpol_pts_ll);
}


// Create polygons for cell-centred data. Here, the polygons are mesh cells
ConservativeSphericalPolygonInterpolation::PolygonArray
ConservativeSphericalPolygonInterpolation::
get_polygons_celldata(FunctionSpace fs, Data::PolygonsData& md) {
    //std::vector<idx_t>& csp2node, std::vector<std::vector<idx_t>>& node2csp,
    //gidx_t& csp_size, std::vector<idx_t>& csp_cell_index, std::vector<idx_t>& csp_index) {
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation: get_polygons_celldata");
    PolygonArray cspolygons(md.csp_size);
    auto mesh = extract_mesh(fs);
    for(idx_t i = 0; i < md.csp_size; ++i) {
        cspolygons[i] = get_csp_celldata(i, mesh, md);
    }
    return cspolygons;
}


// Create polygons for cell-vertex data. Here, the polygons are subcells of mesh cells created as
// 	 (cell_centre, edge_centre, cell_vertex, edge_centre)
// additionally, subcell-to-node and node-to-subcells mapping are computed
ConservativeSphericalPolygonInterpolation::PolygonArray
ConservativeSphericalPolygonInterpolation::
get_polygons_nodedata(FunctionSpace fs, Data::PolygonsData& md) {
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation: get_polygons_nodedata");
    PolygonArray cspolygons(md.csp_size);
    auto mesh = extract_mesh(fs);
    for(idx_t csp_id = 0; csp_id < md.csp_size; ++csp_id) {
        cspolygons[csp_id] = get_csp_nodedata(csp_id, mesh, md);
    }
    return cspolygons;
}


void ConservativeSphericalPolygonInterpolation::do_setup_impl(const Grid& src_grid, const Grid& tgt_grid) {
    ATLAS_TRACE("ConservativeMethod: do_setup_impl");
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
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation: do_setup with grids");

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
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation: do_setup with function spaces");

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
            Log::debug() << "Could not find matrix in cache" << std::endl;
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
    ATLAS_TRACE("ConservativeSphericalPolygonInterpolation: do_setup with function spaces");
    ATLAS_ASSERT(src_fs);
    ATLAS_ASSERT(tgt_fs);

    bool compute_cache = data_->src_.points.empty();

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

    auto& src = sharable_data_->src_;
    auto& tgt = sharable_data_->tgt_;
    auto& timings = sharable_data_->timings;

    if (mpi::size() > 1) {
        // we need src_halo_size >= 2, tgt_halo_size >= 0 for CellColumns
        // if target is NodeColumns, we need:
        //      tgt_halo_size >= 1 and
        //      src_halo_size large enough to cover the the target halo cells in the first row
        int halo_size = 0;
        src_mesh_.metadata().get("halo", halo_size);
        if (halo_size < 2 and order_ == 2) {
            Log::warning() << "The halo size on source mesh should be at least 2 for the 2nd order CSP interpolation." << std::endl;
        }
        if (not tgt_cell_data_) {
            Log::warning() << "The source cells should cover the half of the first target halo row." << std::endl;
        }
        halo_size = 0;
        tgt_mesh_.metadata().get("halo", halo_size);
        if (not tgt_cell_data_ and halo_size == 0) {
            Log::warning() << "The halo size on target mesh should be at least 1 for the target NodeColumns." << std::endl;
        }
    }

    PolygonArray src_csp;
    if (compute_cache) {
        ATLAS_TRACE_SCOPE("Get source polygons") {
            StopWatch stopwatch;
            stopwatch.start();
            init_polygons_data(src_fs_, src);
            src_csp = get_polygons(src_fs_, src);
            stopwatch.stop();
            remap_stat_.memory[Statistics::MEM_SRC_PLG] = memory_of(src_csp);
            timings.source_polygons_assembly = stopwatch.elapsed();
            remap_stat_.counts[Statistics::NUM_SRC_PLG] = src_csp.size();
        }

        ATLAS_TRACE_SCOPE("Get target polygons") {
            StopWatch stopwatch;
            stopwatch.start();
            init_polygons_data(tgt_fs_, tgt);
            stopwatch.stop();
            timings.target_polygons_assembly = stopwatch.elapsed();
            remap_stat_.counts[Statistics::NUM_TGT_PLG] = tgt.csp_size;
        }
    }

    n_spoints_ = src_fs_.size();
    n_tpoints_ = tgt_fs_.size();

    if (compute_cache) {
        intersect_polygons(src_csp);

        // These arrays correspond to functionspace!
        auto& src_points = src.points;
        auto& src_areas  = src.areas;
        src_points.resize(n_spoints_);
        src_areas.resize(n_spoints_);
        {
            ATLAS_TRACE("Store src_areas and src_point");
            if (src.cell_data) {
                for (idx_t scsp_id = 0; scsp_id < src.csp_size; ++scsp_id) {
                    const idx_t spt   = csp_to_cell(scsp_id, src);
                    const auto& s_csp = src_csp[scsp_id];
                    src_points[spt]   = s_csp.centroid();
                    src_areas[spt]    = s_csp.area();
                }
            }
            else {
                auto& src_node2csp = src.node2csp;
                const auto lonlat   = array::make_view<double, 2>(src_mesh_.nodes().lonlat());
                for (idx_t spt = 0; spt < n_spoints_; ++spt) {
                    if (src_node2csp[spt].size() == 0) {
                        // this is a node to which no subpolygon is associated
                        // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                        auto p = PointLonLat{lonlat(spt, 0), lonlat(spt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, src_points[spt]);
                    }
                    else {
                        // .. in the other case, start computing the barycentre
                        src_points[spt] = PointXYZ{0., 0., 0.};
                    }
                    src_areas[spt] = 0.;
                    for (idx_t isubcell = 0; isubcell < src_node2csp[spt].size(); ++isubcell) {
                        idx_t subcell     = src_node2csp[spt][isubcell];
                        const auto& s_csp = src_csp[subcell];
                        src_areas[spt] += s_csp.area();
                        src_points[spt] = src_points[spt] + PointXYZ::mul(s_csp.centroid(), s_csp.area());
                    }
                    double src_point_norm = PointXYZ::norm(src_points[spt]);
                    if (src_point_norm == 0.) {
                        // Probably encountered degenerate subpolygons that lead to zero area. Can happen sometimes in ORCA grid
                        auto p = PointLonLat{lonlat(spt, 0), lonlat(spt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, src_points[spt]);
                        src_point_norm = PointXYZ::norm(src_points[spt]);
                    }
                    src_points[spt] = PointXYZ::div(src_points[spt], src_point_norm);
                }
            }
        }
        auto& tgt_points = tgt.points;
        auto& tgt_areas  = tgt.areas;
        tgt_points.resize(n_tpoints_);
        tgt_areas.resize(n_tpoints_);
        {
            ATLAS_TRACE("Store src_areas and src_point");
            if (tgt_cell_data_) {
                for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
                    auto tgt_csp = get_csp_celldata(tpt, tgt_mesh_, tgt);
                    const auto& t_csp = tgt_csp;
                    tgt_points[tpt]  = t_csp.centroid();
                    tgt_areas[tpt]  = t_csp.area();
                }
            }
            else {
                auto& tgt_node2csp = tgt.node2csp;
                const auto lonlat   = array::make_view<double, 2>(tgt_mesh_.nodes().lonlat());
                for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
                    if (tgt_node2csp[tpt].size() == 0) {
                        // this is a node to which no subpolygon is associated
                        // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                        auto p = PointLonLat{lonlat(tpt, 0), lonlat(tpt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, tgt_points[tpt]);
                    }
                    else {
                        // .. in the other case, start computing the barycentre
                        tgt_points[tpt] = PointXYZ{0., 0., 0.};
                    }
                    tgt_areas[tpt] = 0.;
                    for (idx_t isubcell = 0; isubcell < tgt_node2csp[tpt].size(); ++isubcell) {
                        idx_t subcell     = tgt_node2csp[tpt][isubcell];
                        auto tgt_csp = get_csp_nodedata(subcell, tgt_mesh_, tgt);
                        const auto& t_csp = tgt_csp;
                        tgt_areas[tpt] += t_csp.area();
                        tgt_points[tpt] = tgt_points[tpt] + PointXYZ::mul(t_csp.centroid(), t_csp.area());
                    }
                    double tgt_point_norm = PointXYZ::norm(tgt_points[tpt]);
                    if (tgt_point_norm == 0.) {
                        // Probably encountered degenerate subpolygons cells with zero area as can happen with ORCA sometimes
                        auto p = PointLonLat{lonlat(tpt, 0), lonlat(tpt, 1)};
                        eckit::geometry::Sphere::convertSphericalToCartesian(1., p, tgt_points[tpt]);
                        tgt_point_norm = PointXYZ::norm(tgt_points[tpt]);
                    }
                    tgt_points[tpt] = PointXYZ::div(tgt_points[tpt], tgt_point_norm);
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
            timings.matrix_assembly = stopwatch.elapsed();
        }
    }

    data_->print(Log::debug());

    if (remap_stat_.intersection) {
        setup_stat();
    }
}


namespace {
    uint64_t xorshift(const uint64_t& n, int i){
    return n^(n>>i);
    }
    uint64_t hash(const uint64_t& n){
    uint64_t p = 0x5555555555555555; // pattern of alternating 0 and 1
    uint64_t c = 17316035218449499591ull;// random uneven integer constant;
    return c*xorshift(p*xorshift(n,32),32);
    }
    uint64_t hash_combine(const uint64_t& seed, const uint64_t& v) {
    uint64_t c = 17316035218449499591ull;// random integer constant;
    return hash(v)^(seed+c);
    }

    using PartitionAndRemoteIndex = std::pair<int, idx_t>;
    struct HashPartitionAndRemoteIndex {
        std::size_t operator()(const std::pair<int, idx_t>& pair) const {
            uint64_t h = 0;
            h = hash_combine(h, pair.first);
            h = hash_combine(h, pair.second);
            return h;
        }
    };
}


void ConservativeSphericalPolygonInterpolation::
build_source_kdtree(util::KDTree<idx_t>& kdt_search, double& max_srccell_rad, const PolygonArray& src_csp) const {
    Log::debug() << "Building KDTree via mesh indices." << std::endl;
    const auto& src = sharable_data_->src_;
    const auto node_part  = array::make_view<int, 1>(src_mesh_.nodes().partition());
    const auto node_ridx  = array::make_indexview<idx_t, 1>(src_mesh_.nodes().remote_index());
    const auto node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto cell_part  = array::make_view<int, 1>(src_mesh_.cells().partition());
    const auto cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto cell_ridx  = array::make_indexview<idx_t, 1>(src_mesh_.cells().remote_index());
    auto mpi_rank = mpi::rank();
    ATLAS_TRACE_SCOPE("build kd-tree for source polygons") {
        std::unordered_set<PartitionAndRemoteIndex, HashPartitionAndRemoteIndex> src_kdtree_set;
        auto consider_src = [&](int part, idx_t ridx, int halo) {
            if (not halo) {
                return true;
            }
            else if (part != mpi_rank) {
                return src_kdtree_set.emplace(part, ridx).second;
            }
            return false;
        };
        if (! src_cell_data_) {
            const auto& node2csp = src.node2csp;
            for (idx_t inode = 0; inode < src_mesh_.nodes().size(); ++inode) {
                const auto& csp_ids = node2csp[inode];
                if (!csp_ids.empty()) {
                    if (consider_src(node_part(inode), node_ridx(inode), node_ghost(inode))) {
                        for (const auto& csp_id : csp_ids) {
                            const auto& s_csp = src_csp[csp_id];
                            kdt_search.insert(s_csp.centroid(), csp_id);
                            max_srccell_rad = std::max(max_srccell_rad, s_csp.radius());
                        }
                    }
                }
            }
        }
        else {
            for (idx_t csp_id = 0; csp_id < src.csp_size; ++csp_id) {
                idx_t icell = csp_to_cell(csp_id, src);
                if (consider_src(cell_part(icell), cell_ridx(icell), cell_halo(icell))) {
                    const auto& s_csp = src_csp[csp_id];
                    kdt_search.insert(s_csp.centroid(), icell);
                    max_srccell_rad = std::max(max_srccell_rad, s_csp.radius());
                }
            }
        }
        kdt_search.build();
    }
}


void ConservativeSphericalPolygonInterpolation::
build_source_kdtree_centroid(util::KDTree<idx_t>& kdt_search, double& max_srccell_rad, const PolygonArray& src_csp) const {
    Log::warning() << "Building KDTree via centroid (obsolete)." << std::endl;

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
    ATLAS_TRACE_SCOPE("build kd-tree for source polygons") {
        for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
            const auto& s_csp = src_csp[scell];
            if (not src_already_in((s_csp.centroid()))) {
                kdt_search.insert(s_csp.centroid(), scell);
                max_srccell_rad = std::max(max_srccell_rad, s_csp.radius());
            }
        }
        kdt_search.build();
    }
}


void ConservativeSphericalPolygonInterpolation::intersect_polygons(const PolygonArray& src_csp) {
    ATLAS_TRACE();
    auto& timings = sharable_data_->timings;
    auto& tgt = sharable_data_->tgt_;

    StopWatch stopwatch;
    util::KDTree<idx_t> kdt_search;
    double max_srccell_rad = 0.;
    {
        stopwatch.start();
        kdt_search.reserve(src_csp.size());
        int build_with_centroids = 0;
        const char* ATLAS_BUILD_KDTREE_CENTROID = ::getenv("ATLAS_BUILD_KDTREE_CENTROID");
        if (ATLAS_BUILD_KDTREE_CENTROID != nullptr) {
                build_with_centroids = std::atof(ATLAS_BUILD_KDTREE_CENTROID);
        }
        if (build_with_centroids) {
            build_source_kdtree_centroid(kdt_search, max_srccell_rad, src_csp);
        }
        else {
            build_source_kdtree(kdt_search, max_srccell_rad, src_csp);
        }
        stopwatch.stop();
    }
    timings.source_kdtree_assembly = stopwatch.elapsed();

    StopWatch stopwatch_kdtree_search;
    StopWatch stopwatch_polygon_intersections;

    enum MeshSizeId
    {
        SRC,
        TGT,
        SRC_TGT_INTERSECT,
        TGT_NONINTERSECT
    };
    enum AreaCoverageId
    {
        TOTAL_TGT,
        MAX_TGT
    };
    std::array<size_t, Statistics::NUM_ENUM_SIZE> num_pol{0, 0, 0, 0, 0};
    std::array<double, 2> area_coverage{0., 0.};

    auto& tgt_iparam_ = sharable_data_->tgt_iparam_;
    auto& src_iparam_ = sharable_data_->src_iparam_;
    auto& tgt_csp_size = tgt.csp_size;
    tgt_iparam_.resize(tgt_csp_size);

    if (validate_ || (limit_ != "none")) {
        src_iparam_.resize(src_csp.size());
    }

    // the worst target polygon coverage for analysis of intersection
    std::pair<idx_t, double> worst_tgt_overcover;
    std::pair<idx_t, double> worst_tgt_undercover;
    worst_tgt_overcover.first = -1;
    worst_tgt_overcover.second = -1.;
    worst_tgt_undercover.first = -1;
    worst_tgt_undercover.second = M_PI; // M_PI is the maximum area of a polygon on the unit sphere

    // NOTE: polygon vertex points at distance < pointsSameEPS will be replaced with one point 
    constexpr double pointsSameEPS = 5.e6 * std::numeric_limits<double>::epsilon();

    std::vector<idx_t> intersection_scsp_ids;
    std::vector<double> intersection_weights;
    std::vector<PointXYZ> intersection_src_centroids;

    bool fpe_for_polygon = Polygon::fpe();
    Polygon::fpe(false);
    bool fpe_disabled = fpe_for_polygon ? atlas::library::disable_floating_point_exception(FE_INVALID) : false;
    auto restore_fpe = [fpe_disabled, fpe_for_polygon] {
        if (fpe_disabled) {
            atlas::library::enable_floating_point_exception(FE_INVALID);
        }
        Polygon::fpe(fpe_for_polygon);
    };

    const auto cell_halo   = array::make_view<const int, 1>(tgt_mesh_.cells().halo());
    const auto node_ghost  = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    auto skip_target = [&] (idx_t tcsp_id) -> bool {
        if (tgt.cell_data) {
            auto tcell = csp_to_cell(tcsp_id, tgt);
            return cell_halo(tcell);
        }
        else {
            auto tnode = tgt.csp2node[tcsp_id];
            return node_ghost(tnode);
        }
    };

    ATLAS_TRACE_SCOPE("intersecting polygons") {
        eckit::Channel blackhole;
        eckit::ProgressTimer progress("Intersecting polygons ", 100, " percent", double(10),
                                    tgt_csp_size > 50 ? Log::info() : blackhole);
        float last_progress_percent = 0.00;
        for (idx_t tcsp_id = 0; tcsp_id < tgt_csp_size; ++tcsp_id) {
            auto tgt_csp = get_tgt_csp(tcsp_id); // This needs to be before 'skip_target' because it fills in tgt.csp2node node-data
                                                 // A bit unfortunate... to be improved!
            if (skip_target(tcsp_id)) {
                continue;
            }
            intersection_scsp_ids.resize(0);
            intersection_weights.resize(0);
            intersection_src_centroids.resize(0);
            const auto& t_csp       = tgt_csp;
            double tgt_cover_area   = 0.;
            if (remap_stat_.timings) {stopwatch_kdtree_search.start(); }
            auto scsp_ids = kdt_search.closestPointsWithinRadius(t_csp.centroid(), t_csp.radius() + max_srccell_rad).payloads();
            if (remap_stat_.timings) {stopwatch_kdtree_search.stop(); }
            for (const auto& scsp_id: scsp_ids) {
                const auto& s_csp = src_csp[scsp_id];
                if (remap_stat_.timings) {stopwatch_polygon_intersections.start(); }
                Polygon csp_i = s_csp.intersect(t_csp, nullptr, pointsSameEPS);
                if (remap_stat_.timings) {stopwatch_polygon_intersections.stop(); }
                double csp_i_area = csp_i.area();
#if ATLAS_BUILD_TYPE_DEBUG
                if (validate_) {
                    int pout;
                    if (inside_vertices(t_csp, s_csp, pout) > 2 && csp_i.area() < 3e-16) {
                        dump_intersection("Zero area intersections with inside_vertices", t_csp, src_csp, scsp_ids);
                    }
                }
#endif
                if (csp_i_area > 0.) {
                    intersection_scsp_ids.emplace_back(scsp_id);
                    intersection_weights.emplace_back(csp_i_area);
                    if (order_ == 2 or not matrix_free_ or not matrixAllocated()) {
                        intersection_src_centroids.emplace_back(csp_i.centroid());
                    }
                    tgt_cover_area += csp_i_area;
                    if (validate_ || (limit_ != "none")) {
                        src_iparam_[scsp_id].csp_ids.emplace_back(tcsp_id);
                        if (validate_) {
                            src_iparam_[scsp_id].weights.emplace_back(csp_i_area);
                        }
#if ATLAS_BUILD_TYPE_DEBUG
                        if (csp_i_area > 1.1 * t_csp.area()) {
                            dump_intersection("Intersection larger than target", t_csp, src_csp, scsp_ids);
                        }
                        if (csp_i_area > 1.1 * s_csp.area()) {
                            dump_intersection("Intersection larger than source", t_csp, src_csp, scsp_ids);
                        }
#endif
                    }
                    if (std::abs(1. - tgt_cover_area / tgt_csp.area()) <= 1e-10) {
                        break;
                    }
                }
            }
            const double tgt_cover_err         = std::abs(t_csp.area() - tgt_cover_area);
            const double tgt_cover_err_percent = 100. * tgt_cover_err / t_csp.area();
            if (validate_ && tgt_cover_err_percent > 0.1) {
                if (mpi::size() == 1) {
                    dump_intersection("Target cell not exactly covered", t_csp, src_csp, scsp_ids);
                    if (remap_stat_.intersection) {
                        area_coverage[TOTAL_TGT] += tgt_cover_err;
                        area_coverage[MAX_TGT] = std::max(area_coverage[MAX_TGT], tgt_cover_err);
                    }
                }
            }
            if (intersection_scsp_ids.size() == 0 and remap_stat_.intersection) {
                num_pol[Statistics::NUM_UNCVR_FULL_TGT]++;
            }
            else if (tgt_cover_err_percent > 0.1 && remap_stat_.intersection) {
                num_pol[Statistics::NUM_UNCVR_PART_TGT]++;
            }
            // if (normalise_intersections_ && tgt_cover_err_percent < 1.) {
            //     double wfactor = t_csp.area() / (tgt_cover_area > 0. ? tgt_cover_area : 1.);
            //     for (idx_t i = 0; i < intersection_weights.size(); i++) {
            //         intersection_weights[i] *= wfactor;
            //     }
            // }
            tgt_iparam_[tcsp_id].csp_ids = intersection_scsp_ids;
            tgt_iparam_[tcsp_id].weights = intersection_weights;
            if (order_ == 2 or not matrix_free_ or not matrixAllocated()) {
                tgt_iparam_[tcsp_id].centroids = intersection_src_centroids;
            }
            if (remap_stat_.intersection) {
                num_pol[Statistics::NUM_INT_PLG] += tgt_iparam_[tcsp_id].weights.size();
            }
            if ( double(tcsp_id) / double(tgt_csp_size) > last_progress_percent ) {
                last_progress_percent += 0.01;
                ++progress;
            }
        }
    }

    restore_fpe();

    timings.polygon_intersections  = stopwatch_polygon_intersections.elapsed();
    timings.source_kdtree_search   = stopwatch_kdtree_search.elapsed();
    num_pol[Statistics::NUM_SRC_PLG]                   = src_csp.size();
    num_pol[Statistics::NUM_TGT_PLG]                   = tgt_csp_size;
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(num_pol.data(), num_pol.size(), eckit::mpi::sum());
    }
    remap_stat_.counts[Statistics::NUM_INT_PLG]   = num_pol[Statistics::NUM_INT_PLG];
    remap_stat_.counts[Statistics::NUM_UNCVR_FULL_TGT] = num_pol[Statistics::NUM_UNCVR_FULL_TGT];
    remap_stat_.counts[Statistics::NUM_UNCVR_PART_TGT] = num_pol[Statistics::NUM_UNCVR_PART_TGT];

    const std::string polygon_intersection_folder = "polygon_intersection/";
    if (validate_ && mpi::rank() == 0) {
        if (filesystem::exists(polygon_intersection_folder)) {
            filesystem::remove_all(polygon_intersection_folder);
        }
        filesystem::create_directory(polygon_intersection_folder);
        Log::warning() << "The worst target polygon coverage on the task 0 is in the folder \e[1mpolygon_intersection\e[0m." << std::endl;
    }

    if (validate_) {
        ATLAS_TRACE_SCOPE("compute errors in coverting target cells with intersections") {
            double geo_err_l1   = 0.;
            double geo_err_linf = -1.;
            for (idx_t tcsp_id = 0; tcsp_id < tgt_csp_size; ++tcsp_id) {
                auto tgt_csp = get_tgt_csp(tcsp_id);
                if (skip_target(tcsp_id)) {
                    continue;
                }
                auto& t_csp = tgt_csp;
                auto& tiparam = tgt_iparam_[tcsp_id];
                double tgt_cover_err = 0.;
                for (idx_t icell = 0; icell < tiparam.weights.size(); ++icell) {
                    tgt_cover_err += tiparam.weights[icell];
                }
                tgt_cover_err = 100. * (tgt_cover_err - t_csp.area()) / t_csp.area();
                if (worst_tgt_overcover.second < tgt_cover_err) {
                    worst_tgt_overcover.second = tgt_cover_err;;
                    worst_tgt_overcover.first = tcsp_id;
                }
                if (worst_tgt_undercover.second > tgt_cover_err) {
                    worst_tgt_undercover.second = tgt_cover_err;
                    worst_tgt_undercover.first = tcsp_id;
                }
                if (std::abs(tgt_cover_err) > 0.5) {
                    std::fstream polygon_intersection_info("polygon_intersection", std::ios::out);
                    PointLonLat centre_ll;
                    eckit::geometry::Sphere::convertCartesianToSpherical(1., t_csp.centroid(), centre_ll);
                    polygon_intersection_info << "WARNING tgt cell " << tcsp_id << " over-covering: \e[1m" << tgt_cover_err << "\e[0m %, cell-centre: "
                        << centre_ll <<"\n";
                    polygon_intersection_info << "source indices: " << tiparam.csp_ids << std::endl;
                    dump_intersection("Target cell not exaclty covered", t_csp, src_csp, tiparam.csp_ids);
                    //dump_polygons_to_json(t_csp, src_csp, tiparam.csp_ids, "bad_polygon", 1.e-16);
                }
                geo_err_l1 += std::abs(0.01 * tgt_cover_err * t_csp.area());
                geo_err_linf = std::max(geo_err_linf, std::abs(0.01 * tgt_cover_err * t_csp.area()));
            }
            ATLAS_TRACE_MPI(ALLREDUCE) {
                mpi::comm().allReduceInPlace(geo_err_l1, eckit::mpi::sum());
                mpi::comm().allReduceInPlace(geo_err_linf, eckit::mpi::max());
            }
            remap_stat_.errors[Statistics::ERR_TGT_INTERSECTPLG_L1]   = geo_err_l1 / unit_sphere_area();
            remap_stat_.errors[Statistics::ERR_TGT_INTERSECTPLG_LINF] = geo_err_linf;
        }
    }

    if (validate_) {
        ATLAS_TRACE_SCOPE("compute and output worst target polygon coverage by its intersections with source polygons") {
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
                auto tcsp_id = worst_tgt_overcover.first;
                auto tgt_csp = get_tgt_csp(tcsp_id);
                dump_polygons_to_json(tgt_csp, pointsSameEPS, src_csp, tgt_iparam_[tcsp_id].csp_ids, polygon_intersection_folder, "worst_target_cell_overcover");
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
                auto tcsp_id = worst_tgt_undercover.first;
                auto tgt_csp = get_tgt_csp(tcsp_id);
                dump_polygons_to_json(tgt_csp, pointsSameEPS, src_csp, tgt_iparam_[tcsp_id].csp_ids, polygon_intersection_folder, "worst_target_cell_undercover");
            }
        }
    }
}


ConservativeSphericalPolygonInterpolation::Triplets ConservativeSphericalPolygonInterpolation::compute_1st_order_triplets() {
    ATLAS_TRACE("ConservativeMethod::setup: build cons-1 interpolant matrix");
    ATLAS_ASSERT(not matrix_free_);
    const auto& tgt_iparam = data_->tgt_iparam_;
    Triplets triplets;
    size_t triplets_size    = 0;
    // determine the size of array of triplets used to define the sparse matrix
    if (tgt_cell_data_) {
        for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
            triplets_size += tgt_iparam[tcsp_id].csp_ids.size();
        }
    }
    else {
        auto& tgt_node2csp = data_->tgt_.node2csp;
        for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
            for (idx_t isubcell = 0; isubcell < tgt_node2csp[tnode].size(); ++isubcell) {
                idx_t subcell = tgt_node2csp[tnode][isubcell];
                triplets_size += tgt_iparam[subcell].weights.size();
            }
        }
    }
    triplets.reserve(triplets_size);

    TargetTriplets target_triplets(normalise_);

    struct ScopedDisableFPE {
        ScopedDisableFPE() {
            fpe_disabled = atlas::library::disable_floating_point_exception(FE_DIVBYZERO);
        }
        ~ScopedDisableFPE() {
            if (fpe_disabled) {
                atlas::library::enable_floating_point_exception(FE_DIVBYZERO);
            }
        }
        bool fpe_disabled;
    } disable_fpe;

    // assemble triplets to define the sparse matrix
    const auto& tgt_areas = data_->tgt_.areas;

    if (tgt_cell_data_ && src_cell_data_) {
        for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
            idx_t tcell = csp_to_cell(tcsp_id, data_->tgt_);
            double inv_tgt_weight = (tgt_areas[tcell] > 0. ? 1. / tgt_areas[tcell] : 0.);
            const auto& iparam = tgt_iparam[tcsp_id];
            target_triplets.reset(tcell);
            for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                idx_t scsp_id = iparam.csp_ids[i_scsp];
                idx_t scell = csp_to_cell(scsp_id, data_->src_);
                double weight = iparam.weights[i_scsp] * inv_tgt_weight;
                target_triplets.add(scell, weight);
            }
            target_triplets.emplace_in(triplets);
        }
    }
    else if (not tgt_cell_data_ && src_cell_data_) {
        auto& tgt_node2csp = data_->tgt_.node2csp;
        for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
            double inv_tgt_weight = (tgt_areas[tnode] > 0. ? 1. / tgt_areas[tnode] : 0.);
            target_triplets.reset(tnode);
            for (idx_t isubcell = 0; isubcell < tgt_node2csp[tnode].size(); ++isubcell) {
                const idx_t subcell = tgt_node2csp[tnode][isubcell];
                const auto& iparam  = tgt_iparam[subcell];
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    idx_t scsp_id = iparam.csp_ids[i_scsp];
                    idx_t scell = csp_to_cell(scsp_id, data_->src_);
                    ATLAS_ASSERT(scsp_id < n_spoints_);
                    volatile double weight = iparam.weights[i_scsp] * inv_tgt_weight;
                    target_triplets.add(scell, weight);
                }
            }
            target_triplets.emplace_in(triplets);
        }
    }
    else if (tgt_cell_data_ && not src_cell_data_) {
        auto& src_csp2node = data_->src_.csp2node;
        for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
            idx_t tcell = csp_to_cell(tcsp_id, data_->tgt_);
            double inv_tgt_weight = (tgt_areas[tcell] > 0. ? 1. / tgt_areas[tcell] : 0.);
            const auto& iparam = tgt_iparam[tcsp_id];
            target_triplets.reset(tcell);
            for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                idx_t scsp_id = iparam.csp_ids[i_scsp];
                idx_t snode   = src_csp2node[scsp_id];
                ATLAS_ASSERT(snode < n_spoints_);
                double weight = iparam.weights[i_scsp] * inv_tgt_weight;
                target_triplets.add(snode, weight);
            }
            target_triplets.emplace_in(triplets);
        }
    }
    else if (not tgt_cell_data_ && not src_cell_data_) {
        auto& tgt_node2csp = data_->tgt_.node2csp;
        auto& src_csp2node = data_->src_.csp2node;
        for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
            double inv_tgt_weight = (tgt_areas[tnode] > 0. ? 1. / tgt_areas[tnode] : 0.);
            target_triplets.reset(tnode);
            for (const auto& tcsp_id: tgt_node2csp[tnode]) {
                const auto& iparam  = tgt_iparam[tcsp_id];
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    idx_t scsp_id = iparam.csp_ids[i_scsp];
                    idx_t snode   = src_csp2node[scsp_id];
                    ATLAS_ASSERT(snode < n_spoints_);
                    double weight = iparam.weights[i_scsp] * inv_tgt_weight;
                    target_triplets.add(snode, weight);
                }
            }
            target_triplets.emplace_in(triplets);
        }
    }
    if (validate_) {
        ATLAS_TRACE("ConservativeMethod::setup: Validate the cons-1 matrix");
        std::vector<double> weight_sum(n_tpoints_);
        for (auto& triplet : triplets) {
            weight_sum[triplet.row()] += triplet.value();
        }
        if (order_ == 1) {
            // first order should not give overshoots
            double eps = 1e4 * std::numeric_limits<double>::epsilon();
            for( auto& triplet : triplets ) {
                if (triplet.value() > 1. + eps or triplet.value() < -eps) {
                    Log::info() << "target point " << triplet.row() << " weight: " << triplet.value() << std::endl;
                }
            }
        }
        auto ghost = array::make_view<int,1>( data_->tgt_.cell_data ? tgt_mesh_.cells().halo() : tgt_mesh_.nodes().ghost() );
        for (size_t row = 0; row < n_tpoints_; ++row) {
            if (ghost(row)) {
                continue;
            }
            if (std::abs(weight_sum[row] - 1.) > 1e-10) {
                Log::info() << "target weight in row " << row << " differs from 1 by " << std::abs(weight_sum[row] - 1.) << std::endl;
            }
        }
    }
    return triplets;
}


ConservativeSphericalPolygonInterpolation::Triplets ConservativeSphericalPolygonInterpolation::compute_2nd_order_triplets() {
    ATLAS_TRACE("ConservativeMethod: build cons-2 interpolant matrix");
    Triplets triplets;
    ATLAS_ASSERT(not matrix_free_);
    const auto& src = data_->src_;
    const auto& tgt = data_->tgt_;
    const auto& src_points = src.points;
    const auto& tgt_iparam_ = data_->tgt_iparam_;
    size_t triplets_size    = 0;
    const auto& tgt_areas = data_->tgt_.areas;

    // struct tripple{idx_t row; idx_t column; double value; };
    // std::vector<tripple> triplets;

    struct ScopedDisableFPE {
        ScopedDisableFPE() {
            fpe_disabled = atlas::library::disable_floating_point_exception(FE_DIVBYZERO);
        }
        ~ScopedDisableFPE() {
            if (fpe_disabled) {
                atlas::library::enable_floating_point_exception(FE_DIVBYZERO);
            }
        }
        bool fpe_disabled;
    } disable_fpe;


    Workspace_get_cell_neighbours w_cell;
    Workspace_get_node_neighbours w_node;
    TargetTriplets target_triplets(normalise_);

    if (tgt_cell_data_) {
        const auto tgt_halo = array::make_view<int, 1>(tgt_mesh_.cells().halo());

        for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
            const auto& iparam = tgt_iparam_[tcsp_id];
            idx_t tcell = csp_to_cell(tcsp_id, tgt);
            if (iparam.csp_ids.size() == 0 || tgt_halo(tcell)) {
                continue;
            }
            if (src_cell_data_) {
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    const idx_t scell = csp_to_cell(iparam.csp_ids[i_scsp], src);
                    const auto src_neighbours = get_cell_neighbours(src_mesh_, scell, w_cell);
                    triplets_size += 2 * src_neighbours.size() + 1;
                }
            }
            else {
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    const idx_t scsp_id = iparam.csp_ids[i_scsp];
                    const idx_t snode = src.csp2node[scsp_id];
                    const auto src_neighbours = get_node_neighbours(src_mesh_, snode, w_node);
                    triplets_size += 2 * src_neighbours.size() + 1; 
                }
            }
        }
        triplets.reserve(triplets_size);
        std::vector<PointXYZ> Rsj;
        std::vector<PointXYZ> Aik;
        for (idx_t tcsp_id = 0; tcsp_id < tgt.csp_size; ++tcsp_id) {
            idx_t tcell = csp_to_cell(tcsp_id, tgt);
            const auto& iparam  = tgt_iparam_[tcsp_id];
            if (iparam.csp_ids.size() == 0 || tgt_halo(tcell)) {
                continue;
            }
            target_triplets.reset(tcell);
            // better conservation after Kritsikis et al. (2017)
            // NOTE: ommited here at cost of conservation due to more involved implementation in parallel
            // TEMP: use barycentre computed based on the source polygon's vertices, rather then
            //       the numerical barycentre based on barycentres of the intersections with target cells.
            // TODO: for a given source cell, collect the centroids of all its intersections with target cells
            //       to compute the numerical barycentre of the cell bases on intersection.
            // PointXYZ Cs = {0., 0., 0.};
            // for (idx_t icell = 0; icell < iparam.csp_ids.size(); ++icell) {
            //     Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
            // }
            // Cs = PointXYZ::normalize(Cs);
            double tcell_area_inv = ( tgt_areas[tcell] > 0.) ? 1. / tgt_areas[tcell] : 0.;
            Aik.resize(iparam.csp_ids.size());
            for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                const PointXYZ& Csk     = iparam.centroids[i_scsp];
                const idx_t scsp_id     = iparam.csp_ids[i_scsp];
                const idx_t spt         = src_cell_data_ ? csp_to_cell(scsp_id, src) : src.csp2node[scsp_id];
                const PointXYZ& Cs      = src_points[spt]; // !!! this is centroid and NOT numerically consistent barycentre with the computed source cell intersections
                Aik[i_scsp]             = Csk - Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk - Cs));
                std::vector<idx_t> src_neighbours;
                if (src_cell_data_) {
                    idx_t scell = spt;
                    src_neighbours = get_cell_neighbours(src_mesh_, scell, w_cell);
                }
                else {
                    idx_t snode    = spt;
                    src_neighbours = get_node_neighbours(src_mesh_, snode, w_node);
                }
                Rsj.resize(src_neighbours.size());
                PointXYZ Rs   = {0., 0., 0.};
                double dual_area_inv = 0.;
                for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                    idx_t sj         = src_neighbours[j];
                    idx_t nsj        = src_neighbours[ next_index(j, src_neighbours.size()) ];
                    const auto& Csj  = src_points[sj];
                    const auto& Cnsj = src_points[nsj];
                    if (Polygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16)) {
                        Rsj[j] = PointXYZ::cross(Cnsj, Csj);
                        dual_area_inv += Polygon({Cs, Csj, Cnsj}).area();
                    }
                    else {
                        Rsj[j] = PointXYZ::cross(Csj, Cnsj);
                        dual_area_inv += Polygon({Cs, Cnsj, Csj}).area();
                    }
                    Rs = Rs + Rsj[j];
                }
                dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.;
                Aik[i_scsp]    = PointXYZ::mul(Aik[i_scsp], iparam.weights[i_scsp] * tcell_area_inv * dual_area_inv);

                if (src_cell_data_) {
                    for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                        idx_t nj  = next_index(j, src_neighbours.size());
                        idx_t sj  = src_neighbours[j];
                        idx_t nsj = src_neighbours[nj];
                        double w =  0.5 * PointXYZ::dot(Rsj[j], Aik[i_scsp]);
                        target_triplets.add(sj, w);
                        target_triplets.add(nsj, w);
                    }
                    idx_t scell = spt;
                    target_triplets.add(scell, iparam.weights[i_scsp] * tcell_area_inv - PointXYZ::dot(Rs, Aik[i_scsp]));
                }
                else {
                    for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                        idx_t nj  = next_index(j, src_neighbours.size());
                        idx_t sj  = src_neighbours[j];
                        idx_t nsj = src_neighbours[nj];
                        double w = 0.5 * PointXYZ::dot(Rsj[j], Aik[i_scsp]);
                        target_triplets.add(sj, w);
                        target_triplets.add(nsj, w);
                    }
                    idx_t snode = spt;
                    target_triplets.add(snode, iparam.weights[i_scsp] * tcell_area_inv - PointXYZ::dot(Rs, Aik[i_scsp]));
                }
            }
            target_triplets.emplace_in(triplets);
        }
    }
    else {  // if ( not tgt_cell_data_ )
        const auto tgt_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
        auto& tgt_node2csp = tgt.node2csp;
        triplets_size = 0;
        for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
            if (tgt_ghost(tnode)) {
                continue;
            }
            for (idx_t i_tcsp = 0; i_tcsp < tgt_node2csp[tnode].size(); ++i_tcsp) {
                const idx_t tcsp_id = tgt_node2csp[tnode][i_tcsp];
                const auto& iparam = tgt_iparam_[tcsp_id];
                if (iparam.csp_ids.size() == 0) {
                    continue;
                }
                if (src_cell_data_) {
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        const idx_t scsp_id = iparam.csp_ids[i_scsp];
                        // const idx_t scell = csp_to_cell(scsp_id, data_->src_);
                        const idx_t scell = scsp_id; // !!! TODO: Does this work with meshes with invalid cells ?
                        const auto src_neighbours = get_cell_neighbours(src_mesh_, scell, w_cell);
                        triplets_size += 2 * src_neighbours.size() + 1;
                    }
                }
                else {
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        const idx_t scsp_id = iparam.csp_ids[i_scsp];
                        const idx_t snode = src.csp2node[scsp_id];
                        const auto src_neighbours = get_node_neighbours(src_mesh_, snode, w_node);
                        triplets_size += 2 * src_neighbours.size() + 1; 
                    }
                }
            }
        }

        triplets.reserve(triplets_size);
        std::vector<PointXYZ> Rsj;
        std::vector<PointXYZ> Aik;
        for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
            if (tgt_ghost(tnode)) {
                continue;
            }
            target_triplets.reset(tnode);
            for (idx_t i_tcsp = 0; i_tcsp < tgt_node2csp[tnode].size(); ++i_tcsp) {
                const idx_t tcsp_id = tgt_node2csp[tnode][i_tcsp];
                const auto& iparam = tgt_iparam_[tcsp_id];
                if (iparam.csp_ids.size() == 0) {
                    continue;
                }
                auto tnode_area_inv = ( tgt_areas[tnode] > 0.) ? 1. / tgt_areas[tnode] : 0.;

                // get the barycentre of the dual cell
                // better conservation
                // PointXYZ Cs = {0., 0., 0.};
                // for ( idx_t isubcell = 0; isubcell < src_.node2csp[snode].size(); ++isubcell ) {
                //     idx_t subcell      = src_.node2csp[snode][isubcell];
                //     const auto& iparam = src_iparam_[subcell];
                //     for ( idx_t icell = 0; icell < iparam.csp_ids.size(); ++icell ) {
                //         Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
                //     }
                //     Cs = PointXYZ::normalize(Cs);
                // }
                Aik.resize(iparam.csp_ids.size());
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    const PointXYZ& Csk     = iparam.centroids[i_scsp];
                    const idx_t scsp_id     = iparam.csp_ids[i_scsp];
                    const idx_t spt         = src_cell_data_ ? csp_to_cell(scsp_id, src) : src.csp2node[scsp_id];
                    const PointXYZ& Cs      = src_points[spt]; // !!! this is NOT barycentre of the subcell in case of NodeColumns
                    Aik[i_scsp]             = Csk - Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk - Cs));
                    double dual_area_inv = 0.;
                    std::vector<idx_t> src_neighbours;
                    if (src_cell_data_) {
                        idx_t scell = spt;
                        src_neighbours = get_cell_neighbours(src_mesh_, scell, w_cell);
                    }
                    else {
                        idx_t snode    = spt;
                        src_neighbours = get_node_neighbours(src_mesh_, snode, w_node);
                    }
                    Rsj.resize(src_neighbours.size());
                    PointXYZ Rs   = {0., 0., 0.};
                    for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                        idx_t sj         = src_neighbours[j];
                        idx_t nsj        = src_neighbours[ next_index(j, src_neighbours.size()) ];
                        const auto& Csj  = src_points[sj];
                        const auto& Cnsj = src_points[nsj];
                        if (Polygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16)) {
                            Rsj[j] = PointXYZ::cross(Cnsj, Csj);
                            dual_area_inv += Polygon({Cs, Csj, Cnsj}).area();
                        }
                        else {
                            Rsj[j] = PointXYZ::cross(Csj, Cnsj);
                            dual_area_inv += Polygon({Cs, Cnsj, Csj}).area();
                        }
                        Rs = Rs + Rsj[j];
                    }
                    dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.;
                    Aik[i_scsp]    = PointXYZ::mul(Aik[i_scsp], iparam.weights[i_scsp] * tnode_area_inv * dual_area_inv);

                    if (src_cell_data_) {
                        for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                            idx_t nj  = next_index(j, src_neighbours.size());
                            idx_t sj  = src_neighbours[j];
                            idx_t nsj = src_neighbours[nj];
                            double w = 0.5 * PointXYZ::dot(Rsj[j], Aik[i_scsp]);
                            target_triplets.add(sj, w);
                            target_triplets.add(nsj, w);
                        }
                        idx_t scell = spt;
                        target_triplets.add(scell, iparam.weights[i_scsp] * tnode_area_inv - PointXYZ::dot(Rs, Aik[i_scsp]));
                    }
                    else {
                        for (idx_t j = 0; j < src_neighbours.size(); ++j) {
                            idx_t nj  = next_index(j, src_neighbours.size());
                            idx_t sj  = src_neighbours[j];
                            idx_t nsj = src_neighbours[nj];
                            double w = 0.5 * PointXYZ::dot(Rsj[j], Aik[i_scsp]);
                            target_triplets.add(sj, w);
                            target_triplets.add(nsj, w);
                        }
                        idx_t snode = spt;
                        target_triplets.add(snode, iparam.weights[i_scsp] * tnode_area_inv - PointXYZ::dot(Rs, Aik[i_scsp]));
                    }
                }
            }
            target_triplets.emplace_in(triplets);
        }
    }
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


PointXYZ ConservativeSphericalPolygonInterpolation::src_gradient_celldata(idx_t scell, const array::ArrayView<double, 1>& src_vals) const{
    const auto& src = data_->src_;
    const auto& src_points = src.points;
    PointXYZ grad            = {0., 0., 0.};
    Workspace_get_cell_neighbours w_cell;
    auto nb_cells = get_cell_neighbours(src_mesh_, scell, w_cell);
    volatile double dual_area_inv     = 0.;
    const PointXYZ& Cs       = src_points[scell];
    for (idx_t j = 0; j < nb_cells.size(); ++j) {
        idx_t nj    = next_index(j, nb_cells.size());
        idx_t sj     = nb_cells[j];
        idx_t nsj    = nb_cells[nj];
        const auto& Csj  = src_points[sj];
        const auto& Cnsj = src_points[nsj];
        double val = 0.5 * (src_vals(nj) + src_vals(nsj)) - src_vals(j);
        bool left_orientation = Polygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16);
        dual_area_inv += (left_orientation ? Polygon({Cs, Csj, Cnsj}).area() : Polygon({Cs, Cnsj, Csj}).area());
        grad = grad + PointXYZ::mul(PointXYZ::cross(Csj, Cnsj), val);
    }
    dual_area_inv = ((dual_area_inv > 0.) ? 1. / dual_area_inv : 0.);
    return PointXYZ::mul(grad, dual_area_inv);
}


PointXYZ ConservativeSphericalPolygonInterpolation::src_gradient_nodedata(idx_t snode, const array::ArrayView<double, 1>& src_vals) const{
    const auto& src = data_->src_;
    const auto& src_points = src.points;
    PointXYZ grad            = {0., 0., 0.};
    Workspace_get_node_neighbours w_node;
    auto nb_nodes = get_node_neighbours(src_mesh_, snode, w_node);
    volatile double dual_area_inv     = 0.;
    const PointXYZ& Cs       = src_points[snode]; // TODO: this is not a good barycentre
    for (idx_t j = 0; j < nb_nodes.size(); ++j) {
        idx_t nj    = next_index(j, nb_nodes.size());
        idx_t sj     = nb_nodes[j];
        idx_t nsj    = nb_nodes[nj];
        const auto& Csj  = src_points[sj];
        const auto& Cnsj = src_points[nsj];
        double val = 0.5 * (src_vals(nj) + src_vals(nsj)) - src_vals(j);
        bool left_orientation = Polygon::GreatCircleSegment(Cs, Csj).inLeftHemisphere(Cnsj, -1e-16);
        dual_area_inv += (left_orientation ? Polygon({Cs, Csj, Cnsj}).area() : Polygon({Cs, Cnsj, Csj}).area());
        grad = grad + PointXYZ::mul(PointXYZ::cross(Csj, Cnsj), val);
    }
    dual_area_inv = ((dual_area_inv > 0.) ? 1. / dual_area_inv : 0.);
    return PointXYZ::mul(grad, dual_area_inv);
}


void ConservativeSphericalPolygonInterpolation::do_execute(const Field& src_field, Field& tgt_field,
                                                           Metadata& metadata) const {
    ATLAS_TRACE("ConservativeMethod: do_execute");
    {
        if (src_field.dirty()) {
            ATLAS_TRACE("halo exchange source");
            src_field.haloExchange();
        }
    }
    StopWatch stopwatch;
    stopwatch.start();

    const auto& tgt_iparam = data_->tgt_iparam_;
    std::vector<PointXYZ> src_grads;
    const auto& tgt_areas = data_->tgt_.areas;
    auto& src_iparam = data_->src_iparam_; // will be modified by the limiter
    const auto& src_points = data_->src_.points;
    const auto src_vals = array::make_view<double, 1>(src_field);
    auto tgt_vals       = array::make_view<double, 1>(tgt_field);

    if (! matrix_free_) {
        Method::do_execute(src_field, tgt_field, metadata);
    }
    else {
        if (order_ == 1) {
            ATLAS_TRACE("matrix_free_order_1");

            // CASE: CELL TO CELL
            if (tgt_cell_data_ && src_cell_data_) {
                for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
                    idx_t tcell = csp_to_cell(tcsp_id, data_->tgt_);
                    double tgt_val = 0.;
                    const auto& iparam = tgt_iparam[tcsp_id];
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        idx_t scsp_id = iparam.csp_ids[i_scsp];
                        idx_t scell   = csp_to_cell(scsp_id, data_->src_);
                        tgt_val += iparam.weights[i_scsp] * src_vals(scell);
                    }
                    if (tgt_areas[tcell] > 0.) {
                        tgt_val /= tgt_areas[tcell];
                    }
                    tgt_vals(tcell) = tgt_val;
                }
            }

            // CASE: NODE TO CELL
            else if (not tgt_cell_data_ && src_cell_data_) {
                auto& tgt_node2csp = data_->tgt_.node2csp;
                for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
                    double tgt_val = 0.;
                    for( const auto& tcsp_id: tgt_node2csp[tnode]) {
                        const auto& iparam  = tgt_iparam[tcsp_id];
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t scell   = csp_to_cell(scsp_id, data_->src_);
                            tgt_val += iparam.weights[i_scsp] * src_vals(scell);
                        }
                    }
                    if (tgt_areas[tnode] > 0.) {
                        tgt_val /= tgt_areas[tnode];
                    }
                    tgt_vals(tnode) = tgt_val;
                }
            }

            // CASE: CELL TO NODE
            else if (tgt_cell_data_ && not src_cell_data_) {
                const auto& src_csp2node = data_->src_.csp2node;
                for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
                    idx_t tcell = csp_to_cell(tcsp_id, data_->tgt_);
                    double tgt_val = 0.;
                    const auto& iparam  = tgt_iparam[tcsp_id];
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        idx_t scsp_id = iparam.csp_ids[i_scsp];
                        idx_t snode   = src_csp2node[scsp_id];
                        tgt_val += iparam.weights[i_scsp] * src_vals(snode);
                    }
                    if (tgt_areas[tcell] > 0.) {
                        tgt_val /= tgt_areas[tcell];
                    }
                    tgt_vals(tcell) = tgt_val;
                }
            }

            // CASE: NODE TO NODE
            else if (not tgt_cell_data_ && not src_cell_data_) {
                const auto& tgt_node2csp = data_->tgt_.node2csp;
                const auto& src_csp2node = data_->src_.csp2node;
                for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
                    double tgt_val = 0.;
                    for( const auto& tcsp_id: tgt_node2csp[tnode]) {
                        const auto& iparam = tgt_iparam[tcsp_id];
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t snode   = src_csp2node[scsp_id];
                            tgt_val += iparam.weights[i_scsp] * src_vals(snode);
                        }
                    }
                    if (tgt_areas[tnode] > 0.) {
                        tgt_val /= tgt_areas[tnode];
                    }
                    tgt_vals(tnode) = tgt_val;
                }
            }
        }
        else if (order_ == 2) {
            ATLAS_TRACE_SCOPE("Compute source gradients") {
                src_grads.resize(src_vals.size());
                if (src_cell_data_) {
                    for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                        src_grads[scell] = src_gradient_celldata(scell, src_vals);
                    }
                }
                else {
                    for (idx_t snode = 0; snode < src_vals.size(); ++snode) {
                        src_grads[snode] = src_gradient_nodedata(snode, src_vals);
                    }
                }
            }
            ATLAS_TRACE("matrix_free_order_2");

            // CASE: CELL TO CELL
            if (tgt_cell_data_ && src_cell_data_){
                for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                    const auto& iparam = tgt_iparam[tcell];
                    double tgt_val = 0.;
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        idx_t scsp_id = iparam.csp_ids[i_scsp];
                        idx_t scell = csp_to_cell(scsp_id, data_->src_);
                        const PointXYZ& src_barycentre = src_points[scell]; // TODO: this is a bad barycentre numerically
                        PointXYZ grad  = src_grads[scell];
                        grad           = grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(grad, src_barycentre));
                        tgt_val += iparam.weights[i_scsp] * (src_vals(scell) + PointXYZ::dot(grad, iparam.centroids[i_scsp] - src_barycentre));
                    }
                    if (tgt_areas[tcell] > 0.) {
                        tgt_val /= tgt_areas[tcell];
                    }
                    tgt_vals(tcell) = tgt_val;
                }
            }

            // CASE: NODE TO CELL
            else if (not tgt_cell_data_ && src_cell_data_) {
                auto& tgt_node2csp = data_->tgt_.node2csp;
                for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
                    double tgt_val = 0.;
                    for( const auto& tcsp_id: tgt_node2csp[tnode]) {
                        const auto& iparam  = tgt_iparam[tcsp_id];
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t scell   = csp_to_cell(scsp_id, data_->src_);
                            const PointXYZ& src_barycentre = src_points[scell]; // TODO: this is a bad barycentre numerically
                            PointXYZ grad  = src_grads[scell];
                            grad           = grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(grad, src_barycentre));
                            tgt_val += iparam.weights[i_scsp] * (src_vals(scell) + PointXYZ::dot(grad, iparam.centroids[i_scsp] - src_barycentre));
                        }
                    }
                    if (tgt_areas[tnode] > 0.) {
                        tgt_val /= tgt_areas[tnode];
                    }
                    tgt_vals(tnode) = tgt_val;
                }
            }

            // CASE: CELL TO NODE
            else if (tgt_cell_data_ && not src_cell_data_) {
                const auto& src_csp2node = data_->src_.csp2node;
                for (idx_t tcsp_id = 0; tcsp_id < data_->tgt_.csp_size; ++tcsp_id) {
                    idx_t tcell = csp_to_cell(tcsp_id, data_->tgt_);
                    double tgt_val = 0.;
                    const auto& iparam  = tgt_iparam[tcsp_id];
                    for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                        idx_t scsp_id = iparam.csp_ids[i_scsp];
                        idx_t snode   = src_csp2node[scsp_id];
                        const PointXYZ& src_barycentre = src_points[snode]; // TODO: this is a bad barycentre numerically
                        PointXYZ grad  = src_grads[snode];
                        grad           = grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(grad, src_barycentre));
                        tgt_val += iparam.weights[i_scsp] * (src_vals(snode) + PointXYZ::dot(grad, iparam.centroids[i_scsp] - src_barycentre));
                        // tgt_val += iparam.weights[i_scsp] * src_vals(snode);
                    }
                    if (tgt_areas[tcell] > 0.) {
                        tgt_val /= tgt_areas[tcell];
                    }
                    tgt_vals(tcell) = tgt_val;
                }
            }

            // CASE: NODE TO NODE
            else if (not tgt_cell_data_ && not src_cell_data_) {
                const auto& tgt_node2csp = data_->tgt_.node2csp;
                const auto& src_csp2node = data_->src_.csp2node;
                for (idx_t tnode = 0; tnode < n_tpoints_; ++tnode) {
                    double tgt_val = 0.;
                    for( const auto& tcsp_id: tgt_node2csp[tnode]) {
                        const auto& iparam = tgt_iparam[tcsp_id];
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t snode   = src_csp2node[scsp_id];
                            const PointXYZ& src_barycentre = src_points[snode]; // TODO: this is a bad barycentre numerically
                            PointXYZ grad  = src_grads[snode];
                            grad           = grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(grad, src_barycentre));
                            tgt_val += iparam.weights[i_scsp] * (src_vals(snode) + PointXYZ::dot(grad, iparam.centroids[i_scsp] - src_barycentre));
                            // tgt_val += iparam.weights[i_scsp] * src_vals(snode);
                        }
                    }
                    if (tgt_areas[tnode] > 0.) {
                        tgt_val /= tgt_areas[tnode];
                    }
                    tgt_vals(tnode) = tgt_val;
                }
            }
        }
    }

    struct SrcActed {
        Indices tcells_done;
    };

    if (order_ == 2 && (limit_ != "none")) {
        Field tgt_lim_field = tgt_fs_.createField<double>();
        auto tgt_lim_vals   = array::make_view<double, 1>(tgt_lim_field);
        for (idx_t tcell = 0; tcell < tgt_lim_vals.size(); ++tcell) {
            tgt_lim_vals(tcell) = 0.;
        }

        std::vector<SrcActed> src_acted_tgt;
        src_acted_tgt.resize(src_vals.size());
        double eps = std::numeric_limits<double>::epsilon();

        if (tgt_cell_data_ && src_cell_data_) {
            for (idx_t tcsp = 0; tcsp < data_->tgt_.csp_size; ++tcsp) {
                const auto& iparam = tgt_iparam[tcsp];
                double smax = std::numeric_limits<double>::min();
                double smin = std::numeric_limits<double>::max();
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    idx_t scsp_id = iparam.csp_ids[i_scsp];
                    idx_t scell   = csp_to_cell(scsp_id, data_->src_);
                    smax = std::max(smax, src_vals(scell));
                    smin = std::min(smin, src_vals(scell));
                }
                idx_t tcell = csp_to_cell(tcsp, data_->tgt_);
                double smin_numadj = smin - (1e6 + 1e12 * std::abs(smin)) * eps;
                double smax_numadj = smax + (1e6 + 1e12 * std::abs(smin)) * eps;
                if (tgt_vals(tcell) < smin_numadj || tgt_vals(tcell) > smax_numadj) {
                    if (limit_ == "zeroslope") {
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t scell   = csp_to_cell(scsp_id, data_->src_);
                            const PointXYZ& src_barycentre = src_points[scell];
                            PointXYZ scell_grad  = src_grads[scell];
                            scell_grad           = scell_grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(scell_grad, src_barycentre));
                            auto& siparam = src_iparam[scsp_id];
                            for (idx_t i_tcsp_collateral = 0; i_tcsp_collateral < siparam.csp_ids.size(); ++i_tcsp_collateral) {
                                auto tcsp_collateral = siparam.csp_ids[i_tcsp_collateral];
                                const auto& iparam_collateral = tgt_iparam[tcsp_collateral];
                                auto tcell_collateral = csp_to_cell(tcsp_collateral, data_->tgt_);
                                // find the index of scell entry in iparam_collateral.csp_ids
                                auto scell_it = std::find(iparam_collateral.csp_ids.begin(), iparam_collateral.csp_ids.end(), scell);
                                idx_t scell_idx = scell_it - iparam_collateral.csp_ids.begin();
                                ATLAS_ASSERT(iparam_collateral.csp_ids[scell_idx] == scell);
                                double tgt_lim_val = iparam_collateral.weights[scell_idx] * PointXYZ::dot(scell_grad, iparam_collateral.centroids[scell_idx] - src_barycentre);
                                if (tgt_areas[tcell_collateral] > 0.) {
                                    tgt_lim_val /= tgt_areas[tcell_collateral];
                                }
                                SrcActed& it = src_acted_tgt[scell];
                                if (std::find(it.tcells_done.begin(), it.tcells_done.end(), tcell_collateral) == it.tcells_done.end()) {
                                    it.tcells_done.push_back(tcell_collateral);
                                    tgt_lim_vals(tcell_collateral) -= tgt_lim_val;
                                }
                            }
                        }
                    }
                    else if (limit_ == "clip") {
                        tgt_vals(tcell) = std::max(smin, tgt_vals(tcell));
                        tgt_vals(tcell) = std::min(smax, tgt_vals(tcell));
                    }
                }
            }
            if (limit_ == "zeroslope") {
                for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                    tgt_vals(tcell) += tgt_lim_vals(tcell);
                }
            }
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }
    else if (limit_ == "none") {
        Log::info() << "No remapping limiter applied." << std::endl;
    }

    stopwatch.stop();
    
    if (remap_stat_.conservation) {
        const auto src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
        const auto src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
        const auto src_node_halo  = array::make_view<int, 1>(src_mesh_.nodes().halo());
        const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
        const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
        const auto tgt_node_halo  = array::make_view<int, 1>(tgt_mesh_.nodes().halo());
        const auto& src_areas   = data_->src_.areas;
        const auto& tgt_areas   = data_->tgt_.areas;

        const auto src_vals = array::make_view<double, 1>(src_field);
        const auto tgt_vals = array::make_view<double, 1>(tgt_field);

        double err_remap_cons     = 0.;
        if (src_cell_data_) {
            for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
                if (src_cell_halo(spt)) {
                    continue;
                }
                err_remap_cons += src_vals(spt) * src_areas[spt];
            }
        }
        else {
            for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
                if (src_node_halo(spt) or src_node_ghost(spt)) {
                    continue;
                }
                err_remap_cons += src_vals(spt) * src_areas[spt];
            }
        }
        if (tgt_cell_data_) {
            for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
                if (tgt_cell_halo(tpt)) {
                    continue;
                }
                err_remap_cons -= tgt_vals(tpt) * tgt_areas[tpt];
            }
        }
        else {
            for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
                if (tgt_node_halo(tpt) or tgt_node_ghost(tpt)) {
                    continue;
                }
                err_remap_cons -= tgt_vals(tpt) * tgt_areas[tpt];
            }
        }
        ATLAS_TRACE_MPI(ALLREDUCE) { mpi::comm().allReduceInPlace(&err_remap_cons, 1, eckit::mpi::sum()); }
        remap_stat_.errors[Statistics::ERR_REMAP_CONS] = err_remap_cons / unit_sphere_area();
    }

    if (remap_stat_.intersection) {
        metadata.set("polygons.number_of_src_polygons", remap_stat_.counts[Statistics::NUM_SRC_PLG]);
        metadata.set("polygons.number_of_tgt_polygons", remap_stat_.counts[Statistics::NUM_TGT_PLG]);
        metadata.set("polygons.number_of_intersections", remap_stat_.counts[Statistics::NUM_INT_PLG]);
        metadata.set("polygons.number_of_full_noncovered_tgt_polygons", remap_stat_.counts[Statistics::NUM_UNCVR_FULL_TGT]);
        metadata.set("polygons.number_of_part_noncovered_tgt_polygons", remap_stat_.counts[Statistics::NUM_UNCVR_PART_TGT]);
        if (validate_) {
            metadata.set("errors.intersections_covering_tgt_cells_sum", remap_stat_.errors[Statistics::ERR_TGT_INTERSECTPLG_L1]);
            metadata.set("errors.intersections_covering_tgt_cells_max", remap_stat_.errors[Statistics::ERR_TGT_INTERSECTPLG_LINF]);
        }
    }

    if (remap_stat_.intersection || remap_stat_.conservation) {
        remap_stat_.fillMetadata(metadata);
    }

    auto& timings = data_->timings;
    remap_stat_.time[Statistics::TIME_SRC_PLG] = timings.source_polygons_assembly;
    remap_stat_.time[Statistics::TIME_TGT_PLG] = timings.source_polygons_assembly;
    remap_stat_.time[Statistics::TIME_KDTREE_BUILD] = timings.source_kdtree_assembly;
    remap_stat_.time[Statistics::TIME_MATRIX] = timings.matrix_assembly;
    remap_stat_.time[Statistics::TIME_INTERS] = timings.polygon_intersections;
    if (remap_stat_.timings) {
        remap_stat_.time[Statistics::TIME_KDTREE_SEARCH] = timings.source_kdtree_search;
        metadata.set("max_timings_in_seconds_per_task.source_kdtree_search", timings.source_kdtree_search);
        remap_stat_.time[Statistics::TIME_INTERP] = timings.polygon_intersections;
        metadata.set("max_timings_in_seconds_per_task.polygon_intersections", timings.polygon_intersections);
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(remap_stat_.time.data(), remap_stat_.time.size(), eckit::mpi::max());
    }metadata.set("max_timings_in_seconds_per_task.source_polygons_assembly", timings.source_polygons_assembly);
    metadata.set("max_timings_in_seconds_per_task.target_polygons_assembly", timings.target_polygons_assembly);
    metadata.set("max_timings_in_seconds_per_task.source_kdtree_assembly", timings.source_kdtree_assembly);
    metadata.set("max_timings_in_seconds_per_task.matrix_assembly", timings.matrix_assembly);
    metadata.set("max_timings_in_seconds_per_task.interpolation", stopwatch.elapsed());

    remap_stat_.memory[Statistics::MEM_MATRIX] = (matrix_free_ ? 0 : matrix().footprint());
    remap_stat_.memory[Statistics::MEM_SRC] = memory_of(data_->src_.points);
    remap_stat_.memory[Statistics::MEM_TGT] = memory_of(data_->tgt_.points);
    remap_stat_.memory[Statistics::MEM_SRC_AREAS] = memory_of(data_->src_.areas);
    remap_stat_.memory[Statistics::MEM_TGT_AREAS] = memory_of(data_->tgt_.areas);
    remap_stat_.memory[Statistics::MEM_SRC_CSP2N] = memory_of(data_->src_.csp2node);
    remap_stat_.memory[Statistics::MEM_SRC_N2CSP] = memory_of(data_->src_.node2csp);
    remap_stat_.memory[Statistics::MEM_SRC_CSP2CI] = memory_of(data_->src_.csp_cell_index);
    remap_stat_.memory[Statistics::MEM_SRC_CSP2C] = memory_of(data_->src_.csp_index);
    remap_stat_.memory[Statistics::MEM_TGT_CSP2N] = memory_of(data_->tgt_.csp2node);
    remap_stat_.memory[Statistics::MEM_TGT_N2CSP] = memory_of(data_->tgt_.node2csp);
    remap_stat_.memory[Statistics::MEM_TGT_CSP2CI] = memory_of(data_->tgt_.csp_cell_index);
    remap_stat_.memory[Statistics::MEM_TGT_CSP2C] = memory_of(data_->tgt_.csp_index);
    remap_stat_.memory[Statistics::MEM_IPARAM] = memory_of(data_->tgt_iparam_);
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(remap_stat_.memory.data(), remap_stat_.memory.size(), eckit::mpi::max());
    }
    metadata.set("max_memory_in_bytes_per_task.matrix", remap_stat_.memory[Statistics::MEM_MATRIX]);
    metadata.set("max_memory_in_bytes_per_task.src_points", remap_stat_.memory[Statistics::MEM_SRC]);
    metadata.set("max_memory_in_bytes_per_task.tgt_points", remap_stat_.memory[Statistics::MEM_TGT]);
    metadata.set("max_memory_in_bytes_per_task.src_areas", remap_stat_.memory[Statistics::MEM_SRC_AREAS]);
    metadata.set("max_memory_in_bytes_per_task.tgt_areas", remap_stat_.memory[Statistics::MEM_TGT_AREAS]);
    metadata.set("max_memory_in_bytes_per_task.src_csp2node", remap_stat_.memory[Statistics::MEM_SRC_CSP2N]);
    metadata.set("max_memory_in_bytes_per_task.tgt_csp2node", remap_stat_.memory[Statistics::MEM_TGT_CSP2N]);
    metadata.set("max_memory_in_bytes_per_task.src_node2csp", remap_stat_.memory[Statistics::MEM_SRC_N2CSP]);
    metadata.set("max_memory_in_bytes_per_task.tgt_node2csp", remap_stat_.memory[Statistics::MEM_TGT_N2CSP]);
    metadata.set("max_memory_in_bytes_per_task.tgt_iparam", remap_stat_.memory[Statistics::MEM_IPARAM]);
    metadata.set("max_memory_in_bytes_per_task.src_csp_cell_index", remap_stat_.memory[Statistics::MEM_SRC_CSP2CI]);
    metadata.set("max_memory_in_bytes_per_task.src_csp_index", remap_stat_.memory[Statistics::MEM_SRC_CSP2C]);
    metadata.set("max_memory_in_bytes_per_task.tgt_csp_cell_index", remap_stat_.memory[Statistics::MEM_TGT_CSP2CI]);
    metadata.set("max_memory_in_bytes_per_task.tgt_csp_index", remap_stat_.memory[Statistics::MEM_TGT_CSP2C]);
}


void ConservativeSphericalPolygonInterpolation::print(std::ostream& out) const {
    out << "ConservativeMethod{";
    out << "order:" << order_;
    int halo = 0;
    src_mesh_.metadata().get("halo", halo);
    out << ", source:" << (src_cell_data_ ? "cells(" : "nodes(") << src_mesh_.grid().name() << ",halo=" << halo << ")";
    tgt_mesh_.metadata().get("halo", halo);
    out << ", target:" << (tgt_cell_data_ ? "cells(" : "nodes(") << tgt_mesh_.grid().name() << ",halo=" << halo << ")";
    out << ", normalise:" << normalise_;
    out << ", limit:" << limit_;
    out << ", matrix_free:" << matrix_free_;
    out << ", statistics.intersection:" << remap_stat_.intersection;
    out << ", statistics.conservation:" << remap_stat_.conservation;
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
    if (! remap_stat_.intersection) {
        Log::warning() << "Please enable statistics.intersection." << std::endl;
        return;
    }
    const auto src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto& src_areas   = data_->src_.areas;
    const auto& tgt_areas   = data_->tgt_.areas;
    double geo_create_err     = 0.;
    double src_tgt_sums[2]    = {0., 0.};
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_areas.size(); ++spt) {
            if (not src_cell_halo(spt)) {
                src_tgt_sums[0] += src_areas[spt];
            }
        }
    }
    else {
        for (idx_t spt = 0; spt < src_areas.size(); ++spt) {
            if (not src_node_ghost(spt)) {
                src_tgt_sums[0] += src_areas[spt];
            }
        }
    }
    const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_areas.size(); ++tpt) {
            if (not tgt_cell_halo(tpt)) {
                src_tgt_sums[1] += tgt_areas[tpt];
            }
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_areas.size(); ++tpt) {
            if (not tgt_node_ghost(tpt)) {
                src_tgt_sums[1] += tgt_areas[tpt];
            }
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) { mpi::comm().allReduceInPlace(src_tgt_sums, 2, eckit::mpi::sum()); }

    remap_stat_.src_area_sum = src_tgt_sums[0];
    remap_stat_.tgt_area_sum = src_tgt_sums[1];
    geo_create_err           = std::abs(src_tgt_sums[0] - src_tgt_sums[1]) / unit_sphere_area();
    remap_stat_.errors[Statistics::ERR_SRCTGT_INTERSECTPLG_DIFF] = geo_create_err;
}


void ConservativeSphericalPolygonInterpolation::Statistics::
compute_accuracy(const Interpolation& interpolation, const Field target, std::function<double(const PointLonLat&)> func, Metadata* metadata) {
    auto tgt_vals             = array::make_view<double, 1>(target);
    auto cachable_data_       = ConservativeSphericalPolygonInterpolation::Cache(interpolation).get();
    auto tgt_mesh_            = extract_mesh(cachable_data_->src_fs_);
    auto tgt_cell_data_       = extract_mesh(cachable_data_->tgt_fs_);
    const auto tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    const auto tgt_node_halo  = array::make_view<int, 1>(tgt_mesh_.nodes().halo());
    const auto& tgt_areas     = cachable_data_->tgt_.areas;
    auto& tgt_points          = cachable_data_->tgt_.points;
    double err_remap_l2       = 0.;
    double err_remap_linf     = 0.;
    if (tgt_cell_data_) {
        size_t ncells = std::min<size_t>(tgt_vals.size(), tgt_mesh_.cells().size());
        for (idx_t tpt = 0; tpt < ncells; ++tpt) {
            ATLAS_ASSERT(tpt < tgt_cell_halo.size());
            if (tgt_cell_halo(tpt)) {
                continue;
            }
            auto p = tgt_points[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas[tpt];
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    else {
        size_t nnodes = std::min<size_t>(tgt_vals.size(), tgt_mesh_.nodes().size());
        for (idx_t tpt = 0; tpt < nnodes; ++tpt) {
            if (tgt_node_ghost(tpt) or tgt_node_halo(tpt)) {
                continue;
            }
            auto p = tgt_points[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas[tpt];
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&err_remap_l2, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&err_remap_linf, 1, eckit::mpi::max());
    }
    errors[Statistics::ERR_REMAP_L2]   = std::sqrt(err_remap_l2 / unit_sphere_area());
    errors[Statistics::ERR_REMAP_LINF] = err_remap_linf;
    if (metadata) {
        metadata->set("errors.to_solution_sum", errors[Statistics::ERR_REMAP_L2]);
        metadata->set("errors.to_solution_max", errors[Statistics::ERR_REMAP_LINF]);
    }
}


void ConservativeSphericalPolygonInterpolation::dump_intersection(const std::string msg,
                                                                  const Polygon& plg_1,
                                                                  const PolygonArray& plg_2_array,
                                                                  const std::vector<idx_t>& plg_2_idx_array) const {
#if PRINT_BAD_POLYGONS
    double plg_1_coverage = 0.;
    std::vector<double> int_areas;
    int_areas.resize(plg_2_idx_array.size());
    for (int i = 0; i < plg_2_idx_array.size(); ++i) {
        const auto plg_2_idx = plg_2_idx_array[i];
        const auto& plg_2    = plg_2_array[plg_2_idx].polygon;
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
            const auto& plg_2    = plg_2_array[plg_2_idx].polygon;
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

ConservativeSphericalPolygonInterpolation::Cache::Cache(std::shared_ptr<InterpolationCacheEntry> entry):
    interpolation::Cache(entry), entry_(dynamic_cast<Data*>(entry.get())) {}


ConservativeSphericalPolygonInterpolation::Cache::Cache(const interpolation::Cache& c):
    interpolation::Cache(c, Data::static_type()), entry_{dynamic_cast<const Data*>(c.get(Data::static_type()))} {}


ConservativeSphericalPolygonInterpolation::Cache::Cache(const Interpolation& interpolation):
    Cache(interpolation::Cache(interpolation)) {}


size_t ConservativeSphericalPolygonInterpolation::Data::footprint() const {
    size_t mem_total{0};
    mem_total += memory_of(src_.points);
    mem_total += memory_of(tgt_.points);
    mem_total += memory_of(src_.areas);
    mem_total += memory_of(tgt_.areas);
    mem_total += memory_of(src_.csp2node);
    mem_total += memory_of(tgt_.csp2node);
    mem_total += memory_of(src_.node2csp);
    mem_total += memory_of(tgt_.node2csp);
    mem_total += memory_of(tgt_iparam_);
    // mem_total += memory_of(src_.csp_index); // TODO need to be added
    // mem_total += memory_of(src_.csp_cell_index);
    // mem_total += memory_of(tgt_.csp_index);
    // mem_total += memory_of(src_.csp_cell_index);
    return mem_total;
}


void ConservativeSphericalPolygonInterpolation::Data::print(std::ostream& out) const {
    out << "Memory usage of ConservativeMethod: " << eckit::Bytes(footprint()) << "\n";
    out << "- src_.points   \t" << eckit::Bytes(memory_of(src_.points)) << "\n";
    out << "- tgt_.points   \t" << eckit::Bytes(memory_of(tgt_.points)) << "\n";
    out << "- src_.areas    \t" << eckit::Bytes(memory_of(src_.areas)) << "\n";
    out << "- tgt_.areas    \t" << eckit::Bytes(memory_of(tgt_.areas)) << "\n";
    out << "- src_.csp2node \t" << eckit::Bytes(memory_of(src_.csp2node)) << "\n";
    out << "- tgt_.csp2node \t" << eckit::Bytes(memory_of(tgt_.csp2node)) << "\n";
    out << "- src_.node2csp \t" << eckit::Bytes(memory_of(src_.node2csp)) << "\n";
    out << "- tgt_.node2csp \t" << eckit::Bytes(memory_of(tgt_.node2csp)) << "\n";
    // out << "- src_.csp_index \t" << eckit::Bytes(memory_of(src_.csp_index)) << "\n";
    // out << "- src_csp_cellindex_ \t" << eckit::Bytes(memory_of(src_.csp_cell_index)) << "\n";
    // out << "- tgt_.csp_index \t" << eckit::Bytes(memory_of(tgt_.csp_index)) << "\n";
    // out << "- tgt_csp_cellindex_ \t" << eckit::Bytes(memory_of(tgt_.csp_cell_index)) << "\n";
    out << "- tgt_iparam_   \t" << eckit::Bytes(memory_of(tgt_iparam_)) << "\n";
}


void ConservativeSphericalPolygonInterpolation::Statistics::fillMetadata(Metadata& metadata) {
    // errors
    if (intersection) {
        metadata.set("errors.intersections_covering_tgt_cells_sum", errors[ERR_TGT_INTERSECTPLG_L1]);
        metadata.set("errors.intersections_covering_tgt_cells_max", errors[ERR_TGT_INTERSECTPLG_LINF]);
        metadata.set("errors.sum_src_.areasminus_sum_tgt_areas", errors[ERR_SRCTGT_INTERSECTPLG_DIFF]);
    }
    if (conservation) {
        metadata.set("errors.conservation_error", errors[ERR_REMAP_CONS]);
    }
    if (accuracy) {
        metadata.set("errors.to_solution_sum", errors[ERR_REMAP_L2]);
        metadata.set("errors.to_solution_max", errors[ERR_REMAP_LINF]);
    }
}


ConservativeSphericalPolygonInterpolation::Statistics::Statistics() {
    std::fill(std::begin(counts), std::end(counts), -1);
    std::fill(std::begin(errors), std::end(errors), -1.);
    std::fill(std::begin(memory), std::end(memory), -1);
    std::fill(std::begin(time), std::end(time), -1.);
    intersection = false;
    timings = false;
    accuracy = false;
    conservation = false;
}


ConservativeSphericalPolygonInterpolation::Statistics::Statistics(const Metadata& metadata): Statistics() {
    if (intersection) {
        metadata.get("polygons.number_of_intersections", counts[NUM_INT_PLG]);
        metadata.get("polygons.number_of_full_noncovered_tgt_polygons", counts[NUM_UNCVR_FULL_TGT]);
        metadata.get("polygons.number_of_part_noncovered_tgt_polygons", counts[NUM_UNCVR_PART_TGT]);
    }
    metadata.get("errors.intersections_covering_tgt_cells_sum", errors[ERR_TGT_INTERSECTPLG_L1]);
    metadata.get("errors.intersections_covering_tgt_cells_max", errors[ERR_TGT_INTERSECTPLG_LINF]);
    metadata.get("errors.sum_src_.areasminus_sum_tgt_areas", errors[ERR_SRCTGT_INTERSECTPLG_DIFF]);
    metadata.get("polygons.number_of_src_polygons", counts[NUM_SRC_PLG]);
    metadata.get("polygons.number_of_tgt_polygons", counts[NUM_TGT_PLG]);
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
