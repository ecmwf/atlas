/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/functionspace.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/util/ConvexSphericalPolygon.h"

namespace atlas {
namespace interpolation {
namespace method {

using Indices = std::vector<idx_t>;


class ConservativeSphericalPolygonInterpolation : public Method {
public:
    struct InterpolationParameters {      // one polygon intersection
        Indices csp_ids;                  // target cells used for intersection
        std::vector<PointXYZ> centroids;  // intersection cell centroids
        std::vector<double> weights;      // intersection cell areas
    };

private:
    class Data : public InterpolationCacheEntry {
    public:
        ~Data() override = default;
        size_t footprint() const override;
        static std::string static_type() { return "ConservativeSphericalPolygonInterpolation"; }
        std::string type() const override { return static_type(); }
        void print(std::ostream& out) const;

    private:
        friend class ConservativeSphericalPolygonInterpolation;

        struct PolygonsData {
            enum class Context { SOURCE, TARGET } context;

            // position and effective area of data points
            std::vector<PointXYZ> points;
            std::vector<double>   areas;

            // indexing of subpolygons
            Indices csp2node;
            std::vector<Indices> node2csp;
            gidx_t csp_size;
            Indices csp_cell_index;
            Indices csp_index;
            bool cell_data;
            PolygonsData(PolygonsData::Context ctx) : context(ctx) {}
        } src_{PolygonsData::Context::SOURCE}, tgt_{PolygonsData::Context::TARGET};

        // Timings
        struct Timings {
            double source_polygons_assembly{0};
            double target_polygons_assembly{0};
            double src_already_in{0};
            double source_kdtree_assembly{0};
            double source_kdtree_search{0};
            double source_polygons_filter{0};
            double polygon_intersections{0};
            double matrix_assembly{0};
            double interpolation{0};
        } timings;

        std::vector<InterpolationParameters> tgt_iparam_;

        // Reconstructible if need be
        FunctionSpace src_fs_;
        FunctionSpace tgt_fs_;
    };


    struct TargetTriplets {
        bool normalise_;
        int t;
        std::map<idx_t, double> tpoint_subweights;
        void add (idx_t s, double w) {
            auto pair = tpoint_subweights.emplace(s, w);
            bool inserted = pair.second;
            if (! inserted) {
                auto it = pair.first;
                it->second += w;
            }
        }
        void emplace_in(Triplets& triplets) {
            if (tpoint_subweights.size()) {
                double wfactor = 1.;
                if (normalise_) {
                    volatile double sum_of_weights{0.};
                    for (auto& p : tpoint_subweights) {
                        sum_of_weights += p.second;
                    }
                    ATLAS_ASSERT(sum_of_weights > 0.);
                    wfactor = 1./sum_of_weights;
                }
                for (auto& p : tpoint_subweights) {
                    triplets.emplace_back(t, p.first, p.second * wfactor);
                }
            }
        }
        void reset(int _t) {
            tpoint_subweights.clear();
            t = _t;
        }
        size_t size() const { return tpoint_subweights.size(); }
        TargetTriplets(bool normalise) : normalise_(normalise) {}
    };

public:
    class Cache final : public interpolation::Cache {
    public:
        Cache() = default;
        Cache(const interpolation::Cache& c);
        Cache(const Interpolation&);

        operator bool() const { return entry_; }
        const Data* get() const { return entry_; }

    private:
        friend class ConservativeSphericalPolygonInterpolation;
        Cache(std::shared_ptr<InterpolationCacheEntry> entry);
        const Data* entry_{nullptr};
    };

    struct Statistics {
        enum Counts {
            NUM_SRC_PLG = 0,  // index, number of source polygons
            NUM_TGT_PLG,      // index, number of target polygons
            NUM_INT_PLG,      // index, number of intersection polygons
            NUM_UNCVR_FULL_TGT,    // index, number of completely non covered target polygons
            NUM_UNCVR_PART_TGT,    // index, number of partially non covered target polygons
            NUM_ENUM_SIZE
        };
        enum Errors {
            ERR_TGT_INTERSECTPLG_L1 = 0,      // see above
            ERR_TGT_INTERSECTPLG_LINF,    // see above
            ERR_SRCTGT_INTERSECTPLG_DIFF,    // index, 1/(unit_sphere.area) ( \sum_{scell} scell.area - \sum{tcell} tcell.area )
            ERR_REMAP_CONS,  // index, error in mass conservation
            ERR_REMAP_L2,    // index, error accuracy for given analytical function
            ERR_REMAP_LINF,  // index, like REMAP_L2 but in L_infinity norm
            ERR_ENUM_SIZE
        };
        enum Timings {
            TIME_SRC_PLG = 0,   // index, max time in second per task to build source polygons
            TIME_TGT_PLG,       // index, max time in second per task to build target polygons
            TIME_KDTREE_BUILD,  // index, max time in second per task to build kdtree of source polygons
            TIME_KDTREE_SEARCH, // index, max time in second per task to compute kdtree searches for all target polygons
            TIME_MATRIX,        // index, max time in second per task to assemble the interpolation matrix from the weights
            TIME_INTERS,        // index, max time in second per task to compute intersection polygons
            TIME_INTERP,        // index, max time in second per task to interpolate a source to a target field
            TIME_ENUM_SIZE
        };
        enum Memory {
            MEM_MATRIX = 0,  // index, max memory size per task of the interpolation matrix
            MEM_SRC,         // index, max memory size per task of source point values
            MEM_TGT,
            MEM_SRC_AREAS,   // index, max memory size per task of source-area array
            MEM_TGT_AREAS,
            MEM_SRC_CSP2N,   // index, max memory size per task of source polygon-to-node index-array
            MEM_SRC_N2CSP,   // index, max memory size per task of source node to polygon index-array
            MEM_SRC_CSP2CI,  // index, max memory size per task of source polygon-to-cell-index lookup-array
            MEM_SRC_CSP2C,   // index, max memory size per task of source polygon-to-cell lookup-array
            MEM_SRC_PLG,     // index, max memory size per task of source polygon array
            MEM_TGT_CSP2N,
            MEM_TGT_N2CSP,
            MEM_TGT_CSP2CI,
            MEM_TGT_CSP2C,
            MEM_IPARAM,      // index, max memory size per task of stored intersection parameters
            MEM_ENUM_SIZE
        };
        std::array<int, NUM_ENUM_SIZE> counts;
        std::array<double, ERR_ENUM_SIZE> errors;
        std::array<size_t, MEM_ENUM_SIZE> memory;
        std::array<double, TIME_ENUM_SIZE> time;

        double tgt_area_sum;
        double src_area_sum;
        bool all;
        bool accuracy;
        bool conservation;
        bool intersection;
        bool timings;
        Metadata metadata;

        Statistics();
        Statistics(const Metadata&);

        void compute_accuracy(const Interpolation& interpolation, const Field target, std::function<double(const PointLonLat&)> func, Metadata* metadata = nullptr);
        void fillMetadata(Metadata&);
    };


public:
    ConservativeSphericalPolygonInterpolation(const Config& = util::NoConfig());

    using Method::do_setup;
    void do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs) override;
    void do_setup(const FunctionSpace& source, const FunctionSpace& target, const interpolation::Cache&) override;
    void do_setup(const Grid& src_grid, const Grid& tgt_grid, const interpolation::Cache&) override;

    void do_execute(const Field& src_field, Field& tgt_field, Metadata&) const override;
    void do_execute(const FieldSet& src_fields, FieldSet& tgt_fields, Metadata&) const override;

    void print(std::ostream& out) const override;

    const FunctionSpace& source() const override { return src_fs_; }
    const FunctionSpace& target() const override { return tgt_fs_; }
    Statistics& statistics() { return remap_stat_; }

    inline const PointXYZ& src_points(size_t id) const { return data_->src_.points[id]; }
    inline const PointXYZ& tgt_points(size_t id) const { return data_->tgt_.points[id]; }

    interpolation::Cache createCache() const override;

private:
    using Polygon = util::ConvexSphericalPolygon;
    using PolygonArray = std::vector<util::ConvexSphericalPolygon>;

    void do_setup_impl(const Grid& src_grid, const Grid& tgt_grid);

    void build_source_kdtree(util::KDTree<idx_t>&, double& max_srccell_rad, const PolygonArray&) const;
    void build_source_kdtree_centroid(util::KDTree<idx_t>&, double& max_srccell_rad, const PolygonArray&) const;
    void intersect_polygons(const PolygonArray& src_csp);
    Triplets compute_1st_order_triplets();
    Triplets compute_2nd_order_triplets();
    void dump_intersection(const std::string, const Polygon& plg_1, const PolygonArray& plg_2_array,
                           const Indices& plg_2_idx_array) const;

    struct Workspace_get_cell_neighbours;
    std::vector<idx_t> get_cell_neighbours(Mesh&, idx_t jcell, Workspace_get_cell_neighbours&) const;

    struct Workspace_get_node_neighbours;
    std::vector<idx_t> get_node_neighbours(Mesh&, idx_t jcell, Workspace_get_node_neighbours&) const;

    void init_polygons_data(FunctionSpace fs, Data::PolygonsData& md);


    Polygon get_csp_celldata(idx_t csp_id, const Mesh& mesh, const Data::PolygonsData& md);
    Polygon get_csp_nodedata(idx_t csp_id, const Mesh& mesh, Data::PolygonsData& md);

    idx_t csp_to_cell(idx_t csp_id, const Data::PolygonsData& md) const {
        if (md.cell_data) {
            return md.csp_index[csp_id];
        }
        else {
            auto iterator_upper_bound = std::upper_bound(md.csp_index.begin(), md.csp_index.end(), csp_id);
            idx_t idx = iterator_upper_bound-1 - md.csp_index.begin();
            return md.csp_cell_index[idx];
        }
    }

    std::pair<idx_t, idx_t> csp_to_cell_and_subcell(idx_t csp_id, const Data::PolygonsData& md) const {
        auto iterator_upper_bound = std::upper_bound(md.csp_index.begin(), md.csp_index.end(), csp_id);
        idx_t idx     = iterator_upper_bound-1 - md.csp_index.begin();
        idx_t cell    = md.csp_cell_index[idx];
        idx_t subcell = csp_id - md.csp_index[idx];
        return std::make_pair(cell, subcell);
    }

    // Polygon get_src_csp(idx_t csp_id) {
    //     if (sharable_data_->src_.cell_data) {
    //         return get_csp_celldata(csp_id, src_mesh_, sharable_data_->src_);
    //     }
    //     else {
    //         return get_csp_nodedata(csp_id, src_mesh_, sharable_data_->src_);
    //     }
    // }
    Polygon get_tgt_csp(idx_t csp_id) {
        if (sharable_data_->tgt_.cell_data) {
            return get_csp_celldata(csp_id, tgt_mesh_, sharable_data_->tgt_);
        }
        else {
            return get_csp_nodedata(csp_id, tgt_mesh_, sharable_data_->tgt_);
        }
    }
    PolygonArray get_polygons_celldata(FunctionSpace, Data::PolygonsData&);
    PolygonArray get_polygons_nodedata(FunctionSpace, Data::PolygonsData&);
    PolygonArray get_polygons(FunctionSpace fs, Data::PolygonsData& md) {
        return md.cell_data ? get_polygons_celldata(fs, md) : get_polygons_nodedata(fs, md);
    }


    int next_index(int current_index, int size) const;
    int prev_index(int current_index, int size) const;

    void setup_stat() const;

private:
    bool validate_;
    bool src_cell_data_;
    bool tgt_cell_data_;
    FunctionSpace src_fs_;
    FunctionSpace tgt_fs_;
    mutable Mesh src_mesh_;
    mutable Mesh tgt_mesh_;
    bool normalise_;
    int order_;
    bool matrix_free_;

    mutable Statistics remap_stat_;

    Cache cache_;                          // Storage of cache if any was passed to constructor
    std::shared_ptr<Data> sharable_data_;  // Storage of new data_, only allocated if cache is empty
    const Data* data_;                     // Read-only access to data, pointing either to cache_ or sharable_data_

    // position and effective area of data points
    idx_t n_spoints_;
    idx_t n_tpoints_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
