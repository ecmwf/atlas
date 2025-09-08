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
        Indices cell_idx;      // target cells used for intersection
        std::vector<PointXYZ> centroids;  // intersection cell centroids
        std::vector<double> weights;  // intersection cell areas
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

        // position and effective area of data points
        std::vector<PointXYZ> src_points_;
        std::vector<PointXYZ> tgt_points_;
        std::vector<double> src_areas_;
        std::vector<double> tgt_areas_;

        // indexing of subpolygons
        Indices src_csp2node_;
        Indices tgt_csp2node_;
        std::vector<Indices> src_node2csp_;
        std::vector<Indices> tgt_node2csp_;

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
        enum Counts
        {
            SRC_PLG = 0,  // index, number of source polygons
            TGT_PLG,      // index, number of target polygons
            INT_PLG,      // index, number of intersection polygons
            UNCVR_FULL_TGT,    // index, number of compelte non covered target polygons
            UNCVR_PART_TGT,    // index, number of partially non covered target polygons
            COUNTS_ENUM_SIZE
        };
        enum Errors
        {
            TGT_INTERSECTPLG_L1,      // see above
            TGT_INTERSECTPLG_LINF,    // see above
            SRCTGT_INTERSECTPLG_DIFF,    // index, 1/(unit_sphere.area) ( \sum_{scell} scell.area - \sum{tcell} tcell.area )
            REMAP_CONS,  // index, error in mass conservation
            REMAP_L2,    // index, error accuracy for given analytical function
            REMAP_LINF,  // index, like REMAP_L2 but in L_infinity norm
            ERRORS_ENUM_SIZE
        };
        std::array<int, COUNTS_ENUM_SIZE> counts;
        std::array<double, ERRORS_ENUM_SIZE> errors;

        double tgt_area_sum;
        double src_area_sum;
        bool all;
        bool accuracy;
        bool conservation;
        bool intersection;
        bool timings;
        Metadata metadata;

        void fillMetadata(Metadata&);

        Statistics();
        Statistics(const Metadata&);

        void compute_accuracy(const Interpolation& interpolation, const Field target, std::function<double(const PointLonLat&)> func, Metadata* metadata = nullptr);
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

    inline const PointXYZ& src_points(size_t id) const { return data_->src_points_[id]; }
    inline const PointXYZ& tgt_points(size_t id) const { return data_->tgt_points_[id]; }

    interpolation::Cache createCache() const override;

private:
    using ConvexSphericalPolygon = util::ConvexSphericalPolygon;
    struct MarkedPolygon {
        ConvexSphericalPolygon polygon;
        int halo_type;
    };
    using MarkedPolygonArray           = std::vector<MarkedPolygon>;

    void do_setup_impl(const Grid& src_grid, const Grid& tgt_grid);

    void intersect_polygons(const MarkedPolygonArray& src_csp);
    Triplets compute_1st_order_triplets();
    Triplets compute_2nd_order_triplets();
    void dump_intersection(const std::string, const ConvexSphericalPolygon& plg_1, const MarkedPolygonArray& plg_2_array,
                           const Indices& plg_2_idx_array) const;
    template <class TargetCellsIDs>
    void dump_intersection(const std::string, const ConvexSphericalPolygon& plg_1, const MarkedPolygonArray& plg_2_array,
                           const TargetCellsIDs& plg_2_idx_array) const;
    std::vector<idx_t> sort_cell_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> sort_node_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> get_cell_neighbours(Mesh&, idx_t jcell) const;

    struct Workspace;

    std::vector<idx_t> get_node_neighbours(Mesh&, idx_t jcell, Workspace&) const;
    void init_csp_index(bool cell_data, FunctionSpace fs, Indices& csp2node, std::vector<Indices>& node2csp,
        gidx_t& csp_size, Indices& csp_cell_index, Indices& csp_index);
    MarkedPolygon get_csp_celldata(idx_t csp_id, const Mesh mesh, gidx_t& csp_index_size, Indices& csp_index);
    MarkedPolygon get_csp_nodedata(idx_t csp_id, const Mesh mesh, Indices& csp2node, std::vector<Indices>& node2csp,
        gidx_t& csp_index_size, Indices& csp_cell_index, Indices& csp_index);
    MarkedPolygon get_src_csp(idx_t csp_id) {
        if (src_cell_data_) {
            return get_csp_celldata(csp_id, src_mesh_, src_csp_size_, src_csp_index_);
        }
        else {
            return get_csp_nodedata(csp_id, src_mesh_, sharable_data_->src_csp2node_, sharable_data_->src_node2csp_,
                src_csp_size_, src_csp_cell_index_, src_csp_index_);
        }
    }
    MarkedPolygon get_tgt_csp(idx_t csp_id) {
        if (tgt_cell_data_) {
            return get_csp_celldata(csp_id, tgt_mesh_, tgt_csp_size_, tgt_csp_index_);
        }
        else {
            return get_csp_nodedata(csp_id, tgt_mesh_, sharable_data_->tgt_csp2node_, sharable_data_->tgt_node2csp_,
                tgt_csp_size_, tgt_csp_cell_index_, tgt_csp_index_);
        }
    }
    MarkedPolygonArray get_polygons_celldata(FunctionSpace, Indices& csp2node,
                                         std::vector<Indices>& node2csp,
                                         gidx_t& csp_index_size, Indices& csp_cell_index, Indices& csp_index);
    MarkedPolygonArray get_polygons_nodedata(FunctionSpace, Indices& csp2node,
                                         std::vector<Indices>& node2csp,
                                         gidx_t& csp_index_size, Indices& csp_cell_index, Indices& csp_index);

    int next_index(int current_index, int size, int offset = 1) const;
    int prev_index(int current_index, int size, int offset = 1) const;

    void setup_stat() const;

private:
    bool validate_;
    bool src_cell_data_;
    bool tgt_cell_data_;
    FunctionSpace src_fs_;
    FunctionSpace tgt_fs_;
    mutable Mesh src_mesh_;
    mutable Mesh tgt_mesh_;
    int normalise_intersections_;
    int order_;
    bool matrix_free_;
    bool statistics_timings_;
    bool statistics_intersection_;
    bool statistics_conservation_;

    mutable Statistics remap_stat_;

    gidx_t src_csp_size_;
    gidx_t tgt_csp_size_;
    std::vector<idx_t> src_csp_cell_index_;
    std::vector<idx_t> src_csp_index_;
    std::vector<idx_t> tgt_csp_cell_index_;
    std::vector<idx_t> tgt_csp_index_;

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
