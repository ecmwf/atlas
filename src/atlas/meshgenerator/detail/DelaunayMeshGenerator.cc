/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <unordered_map>

#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildConvexHull3D.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/actions/ExtendNodesGlobal.h"
#include "atlas/meshgenerator/detail/DelaunayMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/projection/Projection.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mesh/ElementType.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

DelaunayMeshGenerator::DelaunayMeshGenerator() = default;

DelaunayMeshGenerator::DelaunayMeshGenerator(const eckit::Parametrisation& p) {
    p.get("part",part_=mpi::rank());
    p.get("reshuffle",reshuffle_=true);
    p.get("remove_duplicate_points",remove_duplicate_points_=true);
}

DelaunayMeshGenerator::~DelaunayMeshGenerator() = default;

void DelaunayMeshGenerator::hash(eckit::Hash& h) const {
    h.add("Delaunay");

    // no other settings
}

void DelaunayMeshGenerator::generate(const Grid& grid, const grid::Distribution& dist, Mesh& mesh) const {
 
    auto build_global_mesh = [&](Mesh& mesh) {
        idx_t nb_nodes = grid.size();
        mesh.nodes().resize(nb_nodes);
        auto xy     = array::make_view<double, 2>(mesh.nodes().xy());
        auto lonlat = array::make_view<double, 2>(mesh.nodes().lonlat());
        auto ghost  = array::make_view<int, 1>(mesh.nodes().ghost());
        auto gidx   = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
        auto part   = array::make_view<int, 1>(mesh.nodes().partition());

        size_t jnode{0};
        Projection projection = grid.projection();
        PointLonLat Pll;
        for (PointXY Pxy : grid.xy()) {
            xy(jnode, size_t(XX)) = Pxy.x();
            xy(jnode, size_t(YY)) = Pxy.y();

            Pll                        = projection.lonlat(Pxy);
            lonlat(jnode, size_t(LON)) = Pll.lon();
            lonlat(jnode, size_t(LAT)) = Pll.lat();

            part(jnode)  = dist.partition(jnode);
            ghost(jnode) = part(jnode) != part_;

            gidx(jnode) = jnode + 1;

            ++jnode;
        }
        mesh::actions::BuildXYZField()(mesh);
        mesh::actions::ExtendNodesGlobal()(grid,mesh);  ///< does nothing if global domain
        mesh::actions::BuildConvexHull3D()(mesh);
        
        auto cells_gidx = array::make_view<gidx_t,1>( mesh.cells().global_index() );
        for (idx_t jelem=0; jelem<mesh.cells().size(); ++jelem) {
            cells_gidx(jelem) = jelem + 1;
        }
    };

    if( dist.nb_partitions() == 1 ) {
        build_global_mesh(mesh);
        return;
    }

    Mesh global_mesh;
    build_global_mesh(global_mesh);

    auto extract_mesh_partition = [&](Mesh& global_mesh, Mesh& mesh) {
        auto g_xy     = array::make_view<double, 2>(global_mesh.nodes().xy());
        auto g_lonlat = array::make_view<double, 2>(global_mesh.nodes().lonlat());
        auto g_ghost  = array::make_view<int, 1>(global_mesh.nodes().ghost());
        auto g_gidx   = array::make_view<gidx_t, 1>(global_mesh.nodes().global_index());
        auto g_part   = array::make_view<int,1>(global_mesh.nodes().partition());

        size_t owned_nodes_count = dist.nb_pts()[part_];

        std::vector<gidx_t> owned_nodes;
        owned_nodes.reserve(1.4*owned_nodes_count);
        for (size_t jnode=0; jnode < global_mesh.nodes().size(); ++ jnode) {
            if (g_ghost(jnode) == 0) {
                owned_nodes.emplace_back(jnode);
            }
        }
        auto& g_node_connectivity = global_mesh.cells().node_connectivity();
        std::set<idx_t> ghost_nodes;
        std::vector<idx_t> owned_elements;
        owned_elements.reserve(1.4*owned_nodes_count);
        std::set<idx_t> element_nodes_uncertainty;
        std::set<idx_t> element_uncertainty;
        constexpr idx_t OWNED = -1;
        constexpr idx_t GHOST = -2;
        constexpr idx_t UNCERTAIN = -3;
        constexpr idx_t CERTAIN = -4;

        auto elem_node_partition = [&](idx_t jelem, idx_t jnode) -> int {
            return g_part(g_node_connectivity(jelem,jnode));
        };

        auto get_elem_ownership = [&](idx_t jelem) -> int {
            int p0 = elem_node_partition(jelem,0);
            int p1 = elem_node_partition(jelem,1);
            int p2 = elem_node_partition(jelem,2);
            if (p0 != part_ && p1 != part_ && p2 != part_) {
                return GHOST;
            }
            if ((p0 == p1 || p0 == p2) && p0 == part_) {
                return OWNED;
            }
            else if (p1 == p2 && p1 == part_) {
                return OWNED;
            }
            else if ( p0 == p1 || p0 == p2 || p1 == p2 ) {
                return CERTAIN;
            }
            return UNCERTAIN;
        };

        auto get_elem_part = [&](idx_t jelem) -> int {
            int p0 = elem_node_partition(jelem,0);
            int p1 = elem_node_partition(jelem,1);
            int p2 = elem_node_partition(jelem,2);
            if (p0 == p1 || p0 == p2) {
                return p0;
            }
            else if (p1 == p2) {
                return p1;
            }
            return UNCERTAIN;
        };


        auto collect_element = [&](idx_t jelem){
            owned_elements.emplace_back(jelem);
            for (idx_t j=0; j<3; ++j) {
                if (elem_node_partition(jelem,j) != part_) {
                    ghost_nodes.insert(g_node_connectivity(jelem,j));
                }
            }
        };

        for (idx_t jelem=0; jelem<global_mesh.cells().size(); ++jelem) {
            idx_t elem_ownership = get_elem_ownership(jelem);
            if (elem_ownership == OWNED) {
                collect_element(jelem);
            }
            else if (elem_ownership == UNCERTAIN) {
                // all three are different
                element_nodes_uncertainty.insert(g_node_connectivity(jelem,0));
                element_nodes_uncertainty.insert(g_node_connectivity(jelem,1));
                element_nodes_uncertainty.insert(g_node_connectivity(jelem,2));
                element_uncertainty.insert(jelem);
            }
        }
        // Log::info() << "element_uncertainty" << std::endl;
        // for( auto& jelem: element_uncertainty ) {
        //     Log::info() << jelem << std::endl;
        // }

        if( element_uncertainty.size() ) {
            std::map<idx_t,std::vector<idx_t>> node2element;
            for (idx_t jelem=0; jelem<global_mesh.cells().size(); ++jelem) {
                for (idx_t jj=0; jj<3; ++jj) {
                    idx_t n = g_node_connectivity(jelem,jj);
                    if (element_nodes_uncertainty.find(n) != element_nodes_uncertainty.end()) {
                        auto it = node2element.find(n);
                        if (it == node2element.end()) {
                            node2element[n].emplace_back(jelem);
                        }
                        else {
                            it->second.emplace_back(jelem);
                        }
                    }
                }
            }
            // Log::info() << "node2element" << std::endl;
            // for( auto& pair: node2element ) {
            //     idx_t jnode = pair.first;
            //     auto& elems = pair.second;
            //     // Log::info() << jnode << " : " << elems << std::endl;
            // }

            auto get_elem_edge = [&](idx_t jelem, idx_t jedge) {
                if (jedge == 0) {
                    return std::array<idx_t,2>{
                        g_node_connectivity(jelem,0),
                        g_node_connectivity(jelem,1),
                    };
                }
                else if(jedge == 1) {
                    return std::array<idx_t,2>{
                        g_node_connectivity(jelem,1),
                        g_node_connectivity(jelem,2),
                    };
                }
                else if(jedge == 2) {
                    return std::array<idx_t,2>{
                        g_node_connectivity(jelem,2),
                        g_node_connectivity(jelem,0),
                    };
                }
                return std::array<idx_t,2>{-1,-1};
            };

            auto get_elem_neighbours = [&](idx_t jelem) -> std::array<idx_t,3> {
                std::array<idx_t,3> elem_neighbours{-1,-1,-1};
                idx_t jneighbour=0;
                for (idx_t jedge=0; jedge<3; ++jedge) {
                    auto edge = get_elem_edge(jelem,jedge);
                    // Log::info() << "jelem,jedge " << jelem << "," << jedge << " : " << edge << "   p: " << g_part(edge[0]) << " " <<   g_part(edge[1]) << std::endl;
                    auto& elem_candidates = node2element.at(edge[0]);
                    for (auto& ielem : elem_candidates) {
                        for (idx_t iedge=0; iedge<3; ++iedge) {
                            auto candidate_edge = get_elem_edge(ielem,iedge);
                            if ( edge[0] == candidate_edge[1] && edge[1] == candidate_edge[0] ) {
                                elem_neighbours[jneighbour++] = ielem;
                                goto next_neighbour;
                            }
                        }
                    }
                    next_neighbour:;
                }
                return elem_neighbours;
            };
            
            for( idx_t jelem : element_uncertainty ) {
                auto elem_neighbours = get_elem_neighbours(jelem);
                idx_t e0 = elem_neighbours[0] >= 0 ? get_elem_part(elem_neighbours[0]) : UNCERTAIN;
                idx_t e1 = elem_neighbours[1] >= 0 ? get_elem_part(elem_neighbours[1]) : UNCERTAIN;
                idx_t e2 = elem_neighbours[2] >= 0 ? get_elem_part(elem_neighbours[2]) : UNCERTAIN;

                idx_t elem_part = UNCERTAIN;
                if (e0 == e1 || e0 == e2) {
                    elem_part = e0;
                }
                else if (e1 == e2) {
                    elem_part = e1;
                }
                else if (e0 != UNCERTAIN) {
                    elem_part = e0;
                }
                else if (e1 != UNCERTAIN) {
                    elem_part = e1;
                }
                else if (e2 != UNCERTAIN) {
                    elem_part = e2;
                }
                if (elem_part == part_) {
                    collect_element(jelem);
                }
            }
        } 

        size_t nb_nodes = owned_nodes.size() + ghost_nodes.size();

        mesh.nodes().resize(nb_nodes);
        auto xy     = array::make_view<double, 2>(mesh.nodes().xy());
        auto lonlat = array::make_view<double, 2>(mesh.nodes().lonlat());
        auto ghost  = array::make_view<int, 1>(mesh.nodes().ghost());
        auto gidx   = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
        auto part   = array::make_view<int, 1>(mesh.nodes().partition());
        auto halo   = array::make_view<int, 1>(mesh.nodes().halo());

        halo.assign(0.);
        std::unordered_map<idx_t,idx_t> from_gnode;
        for (idx_t jnode=0; jnode<owned_nodes.size(); ++jnode) {
            // ATLAS_DEBUG_VAR(jnode);
            idx_t gnode = owned_nodes[jnode];
            from_gnode[gnode] = jnode;
            // ATLAS_DEBUG_VAR(gnode);

            xy(jnode,idx_t(XX)) = g_xy(gnode,idx_t(XX));
            xy(jnode,idx_t(YY)) = g_xy(gnode,idx_t(YY));
            lonlat(jnode,idx_t(XX)) = g_lonlat(gnode,idx_t(XX));
            lonlat(jnode,idx_t(YY)) = g_lonlat(gnode,idx_t(YY));
            ghost(jnode) = 0;
            gidx(jnode) = g_gidx(gnode);
            part(jnode) = part_;
        }
        idx_t jnode = owned_nodes.size();
        for (idx_t gnode: ghost_nodes) {
            ATLAS_ASSERT(jnode < nb_nodes);
            from_gnode[gnode] = jnode;
            xy(jnode,idx_t(XX))     = g_xy(gnode,idx_t(XX));
            xy(jnode,idx_t(YY))     = g_xy(gnode,idx_t(YY));
            lonlat(jnode,idx_t(XX)) = g_lonlat(gnode,idx_t(XX));
            lonlat(jnode,idx_t(YY)) = g_lonlat(gnode,idx_t(YY));
            ghost(jnode) = 1;
            gidx(jnode) = g_gidx(gnode);
            part(jnode) = dist.partition(gnode);
            ++jnode;
        }


        mesh.cells().add(new mesh::temporary::Triangle(), owned_elements.size());
        auto& node_connectivity  = mesh.cells().node_connectivity();
        auto cell_gidx          = array::make_view<gidx_t, 1>(mesh.cells().global_index());
        auto cell_part          = array::make_view<int, 1>(mesh.cells().partition());

        for (idx_t jelem=0; jelem<owned_elements.size(); ++jelem) {
            idx_t gelem = owned_elements[jelem];
            std::array<idx_t,3> triag_nodes {
                from_gnode[g_node_connectivity(gelem,0)],
                from_gnode[g_node_connectivity(gelem,1)],
                from_gnode[g_node_connectivity(gelem,2)]
            };
            node_connectivity.set(jelem, triag_nodes.data());
            cell_gidx(jelem) = gelem+1;
            cell_part(jelem) = part_;
        }
    };

    extract_mesh_partition(global_mesh, mesh);

    setGrid(mesh, grid, dist.type());
}

void DelaunayMeshGenerator::generate(const Grid& g, Mesh& mesh) const {
    generate( g, grid::Distribution{g}, mesh);
}

namespace {
static MeshGeneratorBuilder<DelaunayMeshGenerator> __delaunay("delaunay");
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
