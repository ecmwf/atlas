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

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <type_traits>

#include "eckit/config/Resource.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

using atlas::functionspace::NodeColumns;
using atlas::util::Metadata;
using atlas::util::Topology;
using eckit::PathName;

namespace atlas {
namespace output {
namespace detail {

namespace {

class GmshFile : public std::ofstream {
public:
    GmshFile(const PathName& file_path, std::ios_base::openmode mode, int part = static_cast<int>(mpi::rank())) {
        PathName par_path(file_path);
        int mpi_size = static_cast<int>(mpi::size());
        if (mpi::size() == 1 || part == -1) {
            std::ofstream::open(par_path.localPath(), mode);
        }
        else {
            if (mpi::rank() == 0) {
                PathName par_path(file_path);
                std::ofstream par_file(par_path.localPath(), std::ios_base::out);
                for (int p = 0; p < mpi_size; ++p) {
                    PathName loc_path(file_path);
                    // loc_path = loc_path.baseName(false) + "_p" + to_str(p) + ".msh";
                    loc_path = loc_path.baseName(false) + ".msh.p" + std::to_string(p);
                    par_file << "Merge \"" << loc_path << "\";" << std::endl;
                }
                par_file.close();
            }
            PathName path(file_path);
            // path = path.dirName() + "/" + path.baseName(false) + "_p" +
            // to_str(part) + ".msh";
            path = path.dirName() + "/" + path.baseName(false) + ".msh.p" + std::to_string(part);
            std::ofstream::open(path.localPath(), mode);
        }
    }
};

enum GmshElementTypes
{
    LINE  = 1,
    TRIAG = 2,
    QUAD  = 3,
    POINT = 15
};

// ----------------------------------------------------------------------------
void write_header_ascii(std::ostream& out) {
    out << "$MeshFormat\n";
    out << "4.1 0 " << sizeof(size_t) << "\n";
    out << "$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void write_header_binary(std::ostream& out) {
    out << "$MeshFormat\n";
    out << "4.1 1 " << sizeof(size_t) << "\n";
    int one = 1;
    out.write(reinterpret_cast<const char*>(&one), sizeof(int));
    out << "\n$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------

namespace {  // anonymous

template <typename T>
array::LocalView<const T, 2> make_level_view(const Field& field, int ndata, int jlev) {
    using namespace array;
    if (field.levels()) {
        if (field.variables()) {
            return make_view<const T, 3>(field).slice(Range::to(ndata), jlev, Range::all());
        }
        else {
            return make_view<const T, 2>(field).slice(Range::to(ndata), jlev, Range::dummy());
        }
    }
    else {
        if (field.variables()) {
            return make_view<const T, 2>(field).slice(Range::to(ndata), Range::all());
        }
        else {
            return make_view<const T, 1>(field).slice(Range::to(ndata), Range::dummy());
        }
    }
}

enum class MissingValuePolicy
{
    SKIP,
    FILL,
    NOT_A_NUMBER,
    PRESERVE
};

using MaskedValuePolicy = MissingValuePolicy;

MissingValuePolicy missing_value_policy(const Metadata& gmsh_options) {
    const std::string policy = gmsh_options.get<std::string>("missing_value.policy");
    if (policy == "skip") {
        return MissingValuePolicy::SKIP;
    }
    if (policy == "fill") {
        return MissingValuePolicy::FILL;
    }
    if (policy == "nan") {
        return MissingValuePolicy::NOT_A_NUMBER;
    }
    if (policy == "preserve") {
        return MissingValuePolicy::PRESERVE;
    }
    ATLAS_THROW_EXCEPTION("Unsupported Gmsh missing_value.policy: " << policy);
}

MaskedValuePolicy masked_value_policy(const Metadata& gmsh_options) {
    const std::string policy = gmsh_options.get<std::string>("masked_value.policy");
    if (policy == "skip") {
        return MaskedValuePolicy::SKIP;
    }
    if (policy == "fill") {
        return MaskedValuePolicy::FILL;
    }
    if (policy == "nan") {
        return MaskedValuePolicy::NOT_A_NUMBER;
    }
    if (policy == "preserve") {
        return MaskedValuePolicy::PRESERVE;
    }
    ATLAS_THROW_EXCEPTION("Unsupported Gmsh masked_value.policy: " << policy);
}

template <typename Value>
Value missing_value_fill(const Metadata& gmsh_options, MissingValuePolicy policy) {
    if (policy != MissingValuePolicy::FILL) {
        return Value{0};
    }
    double fill = 0.;
    gmsh_options.get("missing_value.fill", fill);
    return static_cast<Value>(fill);
}

template <typename Value>
Value masked_value_fill(const Metadata& gmsh_options, MaskedValuePolicy policy) {
    if (policy != MaskedValuePolicy::FILL) {
        return Value{0};
    }
    double fill = 0.;
    gmsh_options.get("masked_value.fill", fill);
    return static_cast<Value>(fill);
}

template <typename Value>
class ValueHandling {
public:
    ValueHandling(const Metadata& gmsh_options, const Field& field, const Field& mask = Field()):
        missing_(field),
        policy_(missing_value_policy(gmsh_options)),
        missing_fill_(missing_value_fill<Value>(gmsh_options, policy_)),
        masked_policy_(masked_value_policy(gmsh_options)),
        masked_fill_(masked_value_fill<Value>(gmsh_options, masked_policy_)) {
        if (mask) {
            ATLAS_ASSERT(mask.contiguous());
            mask_ = mask.array().host_data<int>();
        }
    }

    bool include(const array::LocalView<const Value, 2>& data, idx_t n) const {
        if (masked(n)) {
            if (masked_policy_ == MaskedValuePolicy::SKIP) {
                return false;
            }
            if (masked_policy_ != MaskedValuePolicy::PRESERVE) {
                return true;
            }
        }
        if (!missing_ || policy_ != MissingValuePolicy::SKIP) {
            return true;
        }
        for (idx_t v = 0; v < data.shape(1); ++v) {
            if (missing_(data(n, v))) {
                return false;
            }
        }
        return true;
    }

    Value value(Value value, idx_t n) const {
        if (masked(n)) {
            if (masked_policy_ == MaskedValuePolicy::FILL) {
                return masked_fill_;
            }
            if (masked_policy_ == MaskedValuePolicy::NOT_A_NUMBER) {
                if constexpr (std::numeric_limits<Value>::has_quiet_NaN) {
                    return std::numeric_limits<Value>::quiet_NaN();
                }
                ATLAS_THROW_EXCEPTION("Gmsh masked_value.policy 'nan' requires a floating-point field");
            }
        }
        if (!missing_ || !missing_(value)) {
            return value;
        }
        if (policy_ == MissingValuePolicy::FILL) {
            return missing_fill_;
        }
        if (policy_ == MissingValuePolicy::NOT_A_NUMBER) {
            if constexpr (std::numeric_limits<Value>::has_quiet_NaN) {
                return std::numeric_limits<Value>::quiet_NaN();
            }
            ATLAS_THROW_EXCEPTION("Gmsh missing_value.policy 'nan' requires a floating-point field");
        }
        return value;
    }

private:
    bool masked(idx_t n) const {
        return mask_ && mask_[n] == 0;
    }

    field::MissingValue missing_;
    MissingValuePolicy policy_;
    Value missing_fill_;
    MaskedValuePolicy masked_policy_;
    Value masked_fill_;
    const int* mask_{nullptr};
};

class NodeInclusionPolicy {
public:
    explicit NodeInclusionPolicy(const mesh::Nodes& nodes): nodes_(nodes.size()) {
        std::iota(nodes_.begin(), nodes_.end(), idx_t{0});
    }

    const std::vector<idx_t>& included() const { return nodes_; }
    size_t size() const { return nodes_.size(); }
    bool empty() const { return nodes_.empty(); }

private:
    std::vector<idx_t> nodes_;
};

enum class ElementInclusion
{
    SKIP,
    OWNED,
    GHOST
};

class ElementInclusionPolicy {
public:
    ElementInclusionPolicy(const Metadata& gmsh_options, int part, bool include_patch, bool coords_is_lonlat,
                           const array::ArrayView<const double, 2>& coords, const mesh::Elements& elements):
        part_(part),
        include_ghost_(gmsh_options.get<bool>("ghost") && gmsh_options.get<bool>("elements")),
        filter_land_water_(gmsh_options.has("water") || gmsh_options.has("land")),
        include_water_(gmsh_options.getBool("water", false)),
        include_land_(gmsh_options.getBool("land", false)),
        include_patch_(include_patch),
        coords_is_lonlat_(coords_is_lonlat),
        filter_edge_ratio_(eckit::Resource<double>("$ATLAS_GMSH_FILTER_EDGE_RATIO", 0.)),
        coords_(coords),
        elements_(elements),
        nb_nodes_(elements.element_type().name() == "Pentagon" ? 4 : elements.node_connectivity().cols()),
        halo_(elements.view<int, 1>(elements.halo())),
        flags_(elements.view<int, 1>(elements.flags())),
        partition_(elements.view<int, 1>(elements.partition())) {}

    ElementInclusion operator()(idx_t elem) const {
        auto topology = Topology::view(flags_(elem));

        if (filter_land_water_ && !(include_water_ && include_land_) &&
            !((include_water_ && topology.check(Topology::WATER)) ||
              (include_land_ && topology.check(Topology::LAND)))) {
            return ElementInclusion::SKIP;
        }
        if (!include_ghost_ && (topology.check(Topology::GHOST) || halo_(elem))) {
            return ElementInclusion::SKIP;
        }
        if (!include_patch_ && topology.check(Topology::PATCH)) {
            return ElementInclusion::SKIP;
        }
        if (topology.check(Topology::INVALID)) {
            return ElementInclusion::SKIP;
        }
        if (coords_is_lonlat_ && nb_nodes_ == 3 && !include_triangle(elements_.node_connectivity(), elem)) {
            return ElementInclusion::SKIP;
        }
        const bool ghost = (topology.check(Topology::GHOST) || halo_(elem)) && partition_(elem) != part_;
        return ghost ? ElementInclusion::GHOST : ElementInclusion::OWNED;
    }

    size_t nb_nodes() const { return nb_nodes_; }

private:
    bool include_triangle(const mesh::BlockConnectivity& connectivity, idx_t elem) const {
        auto x0            = coords_(connectivity(elem, 0), LON);
        auto x1            = coords_(connectivity(elem, 1), LON);
        auto x2            = coords_(connectivity(elem, 2), LON);
        auto y0            = coords_(connectivity(elem, 0), LAT);
        auto y1            = coords_(connectivity(elem, 1), LAT);
        auto y2            = coords_(connectivity(elem, 2), LAT);
        auto triangle_area = (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1)) * 0.5;
        if (triangle_area <= 0.) {
            return false;
        }
        if (filter_edge_ratio_ <= 0.) {
            return true;
        }
        auto d10  = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
        auto d21  = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
        auto d02  = (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2);
        auto dmin = std::min(d10, std::min(d21, d02));
        auto dmax = std::max(d10, std::max(d21, d02));
        return dmax <= filter_edge_ratio_ * dmin;
    }

    int part_;
    bool include_ghost_;
    bool filter_land_water_;
    bool include_water_;
    bool include_land_;
    bool include_patch_;
    bool coords_is_lonlat_;
    double filter_edge_ratio_;
    array::ArrayView<const double, 2> coords_;
    const mesh::Elements& elements_;
    size_t nb_nodes_;
    array::LocalView<const int, 1> halo_;
    array::LocalView<const int, 1> flags_;
    array::LocalView<const int, 1> partition_;
};

int field_vars(int nvars) {
    if (nvars == 1) {
        return 1;
    }
    if (nvars <= 3) {
        return 3;
    }
    if (nvars == 4 || nvars == 9) {
        return 9;
    }
    ATLAS_NOTIMPLEMENTED;
}

template <typename Value, typename ValueHandling>
typename std::remove_const<Value>::type output_value(const array::LocalView<Value, 2>& data,
                                                     const ValueHandling& value_handling, idx_t n, int component) {
    using value_type     = typename std::remove_const<Value>::type;
    const int input_vars = data.shape(1);
    int input_component  = component;
    if (input_vars == 4) {
        const int row = component / 3;
        const int col = component % 3;
        if (row >= 2 || col >= 2) {
            return value_type{0};
        }
        input_component = row * 2 + col;
    }
    if (input_component >= input_vars) {
        return value_type{0};
    }
    return value_handling.value(data(n, input_component), n);
}

template <typename Value>
void write_ascii_value(std::ostream& out, Value value) {
    if constexpr (std::is_floating_point_v<Value>) {
        if (std::isnan(value)) {
            out << "nan";
            return;
        }
    }
    out << value;
}

template <typename Value, typename GlobalIndex, typename ValueHandling>
void write_level(std::ostream& out, GlobalIndex gidx, const array::LocalView<Value, 2>& data,
                 const ValueHandling& value_handling) {
    const int nvars = field_vars(data.shape(1));
    for (idx_t n = 0; n < data.shape(0); ++n) {
        if (!value_handling.include(data, n)) {
            continue;
        }
        out << gidx(n);
        for (int v = 0; v < nvars; ++v) {
            out << " ";
            write_ascii_value(out, output_value(data, value_handling, n, v));
        }
        out << "\n";
    }
}

template <typename Value, typename GlobalIndex, typename ValueHandling>
void write_level_binary(std::ostream& out, GlobalIndex gidx, const array::LocalView<Value, 2>& data,
                        const ValueHandling& value_handling) {
    const int nvars = field_vars(data.shape(1));
    for (idx_t n = 0; n < data.shape(0); ++n) {
        if (!value_handling.include(data, n)) {
            continue;
        }
        const int tag = static_cast<int>(gidx(n));
        out.write(reinterpret_cast<const char*>(&tag), sizeof(tag));
        for (int v = 0; v < nvars; ++v) {
            const double value = static_cast<double>(output_value(data, value_handling, n, v));
            out.write(reinterpret_cast<const char*>(&value), sizeof(value));
        }
    }
}

template <typename Value, typename ValueHandling>
idx_t count_output_values(const array::LocalView<const Value, 2>& data, const ValueHandling& value_handling) {
    idx_t count = 0;
    for (idx_t n = 0; n < data.shape(0); ++n) {
        count += value_handling.include(data, n);
    }
    return count;
}

std::vector<int> get_levels(int nlev, const Metadata& gmsh_options) {
    std::vector<int> lev;
    std::vector<int> gmsh_levels;
    gmsh_options.get("levels", gmsh_levels);
    if (gmsh_levels.empty() || nlev == 1) {
        lev.resize(nlev);
        for (int ilev = 0; ilev < nlev; ++ilev) {
            lev[ilev] = ilev;
        }
    }
    else {
        lev = gmsh_levels;
    }
    return lev;
}

std::string field_lev(const Field& field, int jlev) {
    if (field.levels()) {
        char str[6]  = {0, 0, 0, 0, 0, 0};
        auto str_len = std::snprintf(str, sizeof(str), "[%03d]", jlev);
        ATLAS_ASSERT(str_len == 5);
        return std::string(str);
    }
    else {
        return std::string();
    }
}

thread_local std::string field_name_prefix;
thread_local Field element_global_index_override;

class FieldNamePrefixScope {
public:
    explicit FieldNamePrefixScope(std::string prefix): previous_(std::move(field_name_prefix)) {
        field_name_prefix = std::move(prefix);
    }

    FieldNamePrefixScope(const FieldNamePrefixScope&)            = delete;
    FieldNamePrefixScope& operator=(const FieldNamePrefixScope&) = delete;

    ~FieldNamePrefixScope() { field_name_prefix = std::move(previous_); }

private:
    std::string previous_;
};

class ElementGlobalIndexScope {
public:
    explicit ElementGlobalIndexScope(const Field& global_index): previous_(element_global_index_override) {
        element_global_index_override = global_index;
    }

    ElementGlobalIndexScope(const ElementGlobalIndexScope&)            = delete;
    ElementGlobalIndexScope& operator=(const ElementGlobalIndexScope&) = delete;

    ~ElementGlobalIndexScope() { element_global_index_override = previous_; }

private:
    Field previous_;
};

std::string output_field_name(const Field& field) {
    return field_name_prefix + field.name();
}

double field_time(const Field& field) {
    return field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
}

int field_step(const Field& field) {
    return field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0;
}

template <typename FunctionSpace>
Field output_mask(const FunctionSpace& function_space, bool gather) {
    if (!function_space.hasMask()) {
        return Field();
    }
    Field mask = function_space.mask();
    if (gather) {
        Field global_mask = function_space.createField(mask, option::global());
        function_space.gather(mask, global_mask);
        return global_mask;
    }
    return mask;
}

Field gmsh_edge_tags(const Mesh& mesh, gidx_t offset) {
    const idx_t nb_edges  = mesh.edges().size();
    auto edge_global_index = array::make_view<gidx_t, 1>(mesh.edges().global_index());

    std::vector<int> counts(mpi::size());
    std::vector<int> displacements(mpi::size());
    constexpr idx_t root = 0;
    mpi::comm().gather(nb_edges, counts, root);

    idx_t nb_global_edges = 0;
    if (mpi::rank() == root) {
        for (idx_t rank = 0; rank < mpi::size(); ++rank) {
            displacements[rank] = nb_global_edges;
            nb_global_edges += counts[rank];
        }
    }

    std::vector<gidx_t> global_edge_tags(nb_global_edges);
    mpi::comm().gatherv(edge_global_index.data(), nb_edges, global_edge_tags.data(), counts.data(),
                        displacements.data(), root);

    if (mpi::rank() == root) {
        std::vector<gidx_t> unique_edge_tags = global_edge_tags;
        std::sort(unique_edge_tags.begin(), unique_edge_tags.end());
        unique_edge_tags.erase(std::unique(unique_edge_tags.begin(), unique_edge_tags.end()),
                               unique_edge_tags.end());
        const bool approximately_compact =
            unique_edge_tags.empty() ||
            unique_edge_tags.back() - unique_edge_tags.front() + 1 <= 2 * unique_edge_tags.size();
        for (gidx_t& tag : global_edge_tags) {
            if (approximately_compact) {
                tag = offset + tag - unique_edge_tags.front() + 1;
            }
            else {
                tag = offset + std::distance(unique_edge_tags.begin(),
                                             std::lower_bound(unique_edge_tags.begin(), unique_edge_tags.end(), tag)) +
                      1;
            }
        }
    }

    Field tags("gmsh_element_tags", array::make_datatype<gidx_t>(), {nb_edges});
    auto local_tags = array::make_view<gidx_t, 1>(tags);
    mpi::comm().scatterv(global_edge_tags.data(), counts.data(), displacements.data(), local_tags.data(), nb_edges,
                         root);

    gidx_t local_max_tag = 0;
    for (idx_t edge = 0; edge < nb_edges; ++edge) {
        local_max_tag = std::max(local_max_tag, local_tags(edge));
    }
    gidx_t max_tag = 0;
    mpi::comm().allReduce(local_max_tag, max_tag, eckit::mpi::max());
    if (max_tag > std::numeric_limits<int>::max()) {
        ATLAS_THROW_EXCEPTION("Gmsh element tag exceeds the supported 32-bit ElementData range");
    }
    return tags;
}

}  // namespace

// ----------------------------------------------------------------------------
template <typename Value>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::NodeColumns& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing NodeColumns field " << field.name() << " defined in NodeColumns..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.nb_nodes(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = array::make_view<gidx_t, 1>(function_space.nodes().global_index());
    Field mask  = output_mask(function_space, gather);
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb =
            function_space.createField<gidx_t>(option::name("gidx_glb") | option::levels(false) | option::global());
        function_space.gather(function_space.nodes().global_index(), gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = std::min<idx_t>(function_space.nb_nodes_global(), field_glb.shape(0));
    }
    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];
        if ((gather && mpi::rank() == 0) || !gather) {
            auto data =
                gather ? make_level_view<Value>(field_glb, ndata, jlev) : make_level_view<Value>(field, ndata, jlev);
            ValueHandling<Value> value_handling(gmsh_options, field, mask);
            const idx_t ndata_output = count_output_values(data, value_handling);

            out << "$NodeData\n";
            out << "1\n";
            out << "\"" << output_field_name(field) << field_lev(field, jlev) << "\"\n";
            out << "1\n";
            out << field_time(field) << "\n";
            out << "4\n";
            out << field_step(field) << "\n";
            out << field_vars(nvars) << "\n";
            out << ndata_output << "\n";
            out << mpi::rank() << "\n";
            if (gmsh_options.get<bool>("ascii")) {
                write_level(out, gidx, data, value_handling);
            }
            else {
                write_level_binary(out, gidx, data, value_handling);
                out << "\n";
            }
            out << "$EndNodeData\n";
        }
    }
}

// ----------------------------------------------------------------------------
template <typename Value>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::NoFunctionSpace& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing field " << field.name() << " defined without functionspace..." << std::endl;

    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = field.shape(0);
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = [](idx_t inode) { return inode + 1; };

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];

        auto data = make_level_view<Value>(field, ndata, jlev);
        ValueHandling<Value> value_handling(gmsh_options, field);
        const idx_t ndata_output = count_output_values(data, value_handling);

        out << "$NodeData\n";
        out << "1\n";
        out << "\"" << output_field_name(field) << field_lev(field, jlev) << "\"\n";
        out << "1\n";
        out << field_time(field) << "\n";
        out << "4\n";
        out << field_step(field) << "\n";
        out << field_vars(nvars) << "\n";
        out << ndata_output << "\n";
        out << mpi::rank() << "\n";
        if (gmsh_options.get<bool>("ascii")) {
            write_level(out, gidx, data, value_handling);
        }
        else {
            write_level_binary(out, gidx, data, value_handling);
            out << "\n";
        }
        out << "$EndNodeData\n";
    }
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

void print_field_lev(char field_lev[], size_t size, int jlev) {
    ATLAS_ASSERT(size > 5);
    std::snprintf(field_lev, size, "[%03d]", jlev);
}

/* unused
void print_field_lev( char field_lev[], long jlev ) {
    std::sprintf( field_lev, "[%03ld]", jlev );
}

void print_field_lev( char field_lev[], unsigned long jlev ) {
    std::sprintf( field_lev, "[%03lu]", jlev );
}
**/

// ----------------------------------------------------------------------------
template <typename DATATYPE>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::StructuredColumns& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing StructuredColumns field " << field.name() << "..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.sizeOwned(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = array::make_view<gidx_t, 1>(function_space.global_index());
    Field mask  = output_mask(function_space, gather);
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb =
            function_space.createField(function_space.global_index(), option::name("gidx_glb") | option::global());
        function_space.gather(function_space.global_index(), gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = field_glb.shape(0);
    }

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev          = lev[ilev];
        char field_lev[6] = {0, 0, 0, 0, 0, 0};

        if (field.levels()) {
            print_field_lev(field_lev, sizeof(field_lev), jlev);
        }

        auto data =
            gather ? make_level_view<DATATYPE>(field_glb, ndata, jlev) : make_level_view<DATATYPE>(field, ndata, jlev);
        ValueHandling<DATATYPE> value_handling(gmsh_options, field, mask);
        const idx_t ndata_output = count_output_values(data, value_handling);


        out << "$NodeData\n";
        out << "1\n";
        out << "\"" << output_field_name(field) << field_lev << "\"\n";
        out << "1\n";
        out << field_time(field) << "\n";
        out << "4\n";
        out << field_step(field) << "\n";
        out << field_vars(nvars) << "\n";
        out << ndata_output << "\n";
        out << mpi::rank() << "\n";
        if (gmsh_options.get<bool>("ascii")) {
            write_level(out, gidx, data, value_handling);
        }
        else {
            write_level_binary(out, gidx, data, value_handling);
            out << "\n";
        }
        out << "$EndNodeData\n";
    }
}
// ----------------------------------------------------------------------------

template <typename DATATYPE, typename FunctionSpace>
void write_field_elems(const Metadata& gmsh_options, const FunctionSpace& function_space,
                       const Field& field, std::ostream& out) {
#if 1
    Log::debug() << "writing element field " << field.name() << "..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.size(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    Field global_index = element_global_index_override ? element_global_index_override : function_space.global_index();
    auto gidx          = array::make_view<gidx_t, 1>(global_index);
    Field mask  = output_mask(function_space, gather);
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb = function_space.template createField<gidx_t>(option::name("gidx_glb") | option::levels(false) |
                                                               option::global());
        function_space.gather(global_index, gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = field_glb.shape(0);
    }

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];

        auto data =
            gather ? make_level_view<DATATYPE>(field_glb, ndata, jlev) : make_level_view<DATATYPE>(field, ndata, jlev);
        ValueHandling<DATATYPE> value_handling(gmsh_options, field, mask);
        const idx_t ndata_output = count_output_values(data, value_handling);


        if ((gather && mpi::rank() == 0) || !gather) {
            out << "$ElementData\n";
            out << "1\n";
            out << "\"" << output_field_name(field) << field_lev(field, jlev) << "\"\n";
            out << "1\n";
            out << field_time(field) << "\n";
            out << "4\n";
            out << field_step(field) << "\n";
            out << field_vars(nvars) << "\n";
            out << ndata_output << "\n";
            out << mpi::rank() << "\n";
            if (gmsh_options.get<bool>("ascii")) {
                write_level(out, gidx, data, value_handling);
            }
            else {
                write_level_binary(out, gidx, data, value_handling);
                out << "\n";
            }
            out << "$EndElementData\n";
        }
    }
#endif
}

// ----------------------------------------------------------------------------
#if 0
template< typename DATA_TYPE >
void write_field_elems(const Metadata& gmsh_options, const FunctionSpace& function_space, const Field& field, std::ostream& out)
{
  Log::info() << "writing field " << field.name() << "..." << std::endl;
  bool gather( gmsh_options.get<bool>("gather") );
  bool binary( !gmsh_options.get<bool>("ascii") );
  int nlev = field.metadata().has("nb_levels") ? field.metadata().get<size_t>("nb_levels") : 1;
  int ndata = field.shape(0);
  int nvars = field.shape(1)/nlev;
  array::ArrayView<gidx_t,1    > gidx ( function_space.field( "glb_idx" ) );
  array::ArrayView<DATA_TYPE> data ( field );
  array::ArrayT<DATA_TYPE> field_glb_arr;
  array::ArrayT<gidx_t   > gidx_glb_arr;
  if( gather )
  {
    mpl::GatherScatter& fullgather = function_space.fullgather();
    ndata = fullgather.glb_dof();
    field_glb_arr.resize(ndata,field.shape(1));
    gidx_glb_arr.resize(ndata);
    array::ArrayView<DATA_TYPE> data_glb( field_glb_arr );
    array::ArrayView<gidx_t,1> gidx_glb( gidx_glb_arr );
    fullgather.gather( gidx, gidx_glb );
    fullgather.gather( data, data_glb );
    gidx = array::ArrayView<gidx_t,1>( gidx_glb_arr );
    data = data_glb;
  }

  double time = field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
  size_t step = field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0 ;

  int nnodes = IndexView<int,2>( function_space.field("nodes") ).shape(1);

  for (int jlev=0; jlev<nlev; ++jlev)
  {
    char field_lev[6] = {0, 0, 0, 0, 0, 0};
    if( field.metadata().has("nb_levels") )
      std::sprintf(field_lev, "[%03d]",jlev);

    out << "$ElementNodeData\n";
    out << "1\n";
    out << "\"" << field.name() << field_lev << "\"\n";
    out << "1\n";
    out << time << "\n";
    out << "4\n";
    out << step << "\n";
    if     ( nvars == 1 ) out << nvars << "\n";
    else if( nvars <= 3 ) out << 3     << "\n";
    out << ndata << "\n";
    out << mpi::rank() << "\n";

    if( binary )
    {
      if( nvars == 1)
      {
        double value;
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out.write(reinterpret_cast<const char*>(&gidx(jelem)),sizeof(int));
          out.write(reinterpret_cast<const char*>(&nnodes),sizeof(int));
          for (size_t n=0; n<nnodes; ++n)
          {
            value = data(jelem,jlev);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double));
          }
        }
      }
      else if( nvars <= 3 )
      {
        double value[3] = {0,0,0};
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              value[v] = data(jelem,jlev*nvars+v);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double)*3);
          }
        }
      }
      out <<"\n";
    }
    else
    {
      if( nvars == 1)
      {
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
            out << " " << data(jelem,jlev);
          out <<"\n";
        }
      }
      else if( nvars <= 3 )
      {
        std::vector<DATA_TYPE> data_vec(3,0.);
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              data_vec[v] = data(jelem,jlev*nvars+v);
            for( int v=0; v<3; ++v)
              out << " " << data_vec[v];
          }
          out <<"\n";
        }
      }
    }
    out << "$EndElementNodeData\n";
  }
}
#endif

// ----------------------------------------------------------------------------

}  // end anonymous namespace

// ----------------------------------------------------------------------------
GmshIO::GmshIO() {
    // which field holds the Nodes
    options.set<std::string>("nodes", "xy");

    // Gather fields to one proc before writing
    options.set<bool>("gather", false);

    // Output of ghost nodes / elements
    options.set<bool>("ghost", false);

    // ASCII format (true) or binary (false)
    options.set<bool>("ascii", true);

    // Output of elements
    options.set<bool>("elements", true);

    // Output of edges
    options.set<bool>("edges", true);

    // Use the zero-based element owner as the elementary entity tag
    options.set<bool>("element_partition_as_entity", false);

    // Levels of fields to use
    options.set<std::vector<long>>("levels", std::vector<long>());

    options.set<std::string>("missing_value.policy", "skip");
    options.set<double>("missing_value.fill", 0.);
    options.set<std::string>("masked_value.policy", "skip");
    options.set<double>("masked_value.fill", 0.);
}

GmshIO::~GmshIO() = default;

Mesh GmshIO::read(const PathName& file_path) const {
    Mesh mesh;
    GmshIO::read(file_path, mesh);
    return mesh;
}

namespace {
mesh::ElementType* make_element_type(int type) {
    if (type == QUAD) {
        return mesh::ElementType::create("Quadrilateral");
    }
    if (type == TRIAG) {
        return mesh::ElementType::create("Triangle");
    }
    if (type == LINE) {
        return mesh::ElementType::create("Line");
    }
    throw_Exception("Element type not supported", Here());
}
}  // namespace

void GmshIO::read(const PathName& file_path, Mesh& mesh) const {
    std::ifstream file;
    file.open(file_path.localPath(), std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        throw_CantOpenFile(file_path);
    }

    std::string line;

    while (line != "$MeshFormat") {
        std::getline(file, line);
    }
    double version;
    int binary;
    int size_of_real;
    file >> version >> binary >> size_of_real;

    while (line != "$Nodes") {
        std::getline(file, line);
    }

    // Create nodes
    idx_t nb_nodes;
    file >> nb_nodes;

    mesh.nodes().resize(nb_nodes);

    mesh::Nodes& nodes = mesh.nodes();

    //nodes.add( Field( "xyz", array::make_datatype<double>(), array::make_shape( nb_nodes, 3 ) ) );

    //    array::ArrayView<double, 2> coords  = array::make_view<double, 2>( nodes.field( "xyz" ) );
    array::ArrayView<double, 2> xy      = array::make_view<double, 2>(nodes.xy());
    array::ArrayView<double, 2> lonlat  = array::make_view<double, 2>(nodes.lonlat());
    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());
    array::ArrayView<int, 1> part       = array::make_view<int, 1>(nodes.partition());
    array::ArrayView<int, 1> ghost      = array::make_view<int, 1>(nodes.ghost());

    std::map<int, int> glb_to_loc;
    int g;
    double x, y, z;
    double xyz[3];
    double xmax        = -std::numeric_limits<double>::max();
    double zmax        = -std::numeric_limits<double>::max();
    gidx_t max_glb_idx = 0;
    while (binary && file.peek() == '\n') {
        file.get();
    }
    for (idx_t n = 0; n < nb_nodes; ++n) {
        if (binary) {
            file.read(reinterpret_cast<char*>(&g), sizeof(int));
            file.read(reinterpret_cast<char*>(&xyz), sizeof(double) * 3);
            x = xyz[XX];
            y = xyz[YY];
            z = xyz[ZZ];
        }
        else {
            file >> g >> x >> y >> z;
        }
        glb_idx(n)     = g;
        xy(n, XX)      = x;
        xy(n, YY)      = y;
        lonlat(n, LON) = x;
        lonlat(n, LAT) = y;
        glb_to_loc[g]  = n;
        part(n)        = 0;
        ghost(n)       = 0;
        max_glb_idx    = std::max(max_glb_idx, static_cast<gidx_t>(g));
        xmax           = std::max(x, xmax);
        zmax           = std::max(z, zmax);
    }
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    int nb_elements = 0;

    while (line != "$Elements") {
        std::getline(file, line);
    }

    file >> nb_elements;

    if (binary) {
        while (file.peek() == '\n') {
            file.get();
        }
        int accounted_elems = 0;
        while (accounted_elems < nb_elements) {
            int header[3];
            int data[100];
            file.read(reinterpret_cast<char*>(&header), sizeof(int) * 3);

            int etype     = header[0];
            size_t netype = header[1];
            size_t ntags  = header[2];
            accounted_elems += netype;
            mesh::Elements* elements;
            if (etype == LINE) {
                size_t jtype = mesh.edges().add(make_element_type(etype), netype);
                elements     = &mesh.edges().elements(jtype);
            }
            else {
                size_t jtype = mesh.cells().add(make_element_type(etype), netype);
                elements     = &mesh.cells().elements(jtype);
            }

            size_t nnodes_per_elem            = elements->element_type().nb_nodes();
            mesh::BlockConnectivity& conn     = elements->node_connectivity();
            array::ArrayView<gidx_t, 1> egidx = array::make_view<gidx_t, 1>(elements->global_index());
            array::ArrayView<int, 1> epart    = array::make_view<int, 1>(elements->partition());

            size_t dsize = 1 + ntags + nnodes_per_elem;
            int part;
            for (size_t e = 0; e < netype; ++e) {
                file.read(reinterpret_cast<char*>(&data), sizeof(int) * dsize);
                part     = 0;
                egidx(e) = data[0];
                epart(e) = part;
                for (size_t n = 0; n < nnodes_per_elem; ++n) {
                    conn.set(e, n, glb_to_loc[data[1 + ntags + n]]);
                }
            }
        }
    }
    else {
        // Find out which element types are inside
        int position = file.tellg();
        std::vector<int> nb_etype(20, 0);
        int elements_max_glb_idx(0);
        int etype;
        for (int e = 0; e < nb_elements; ++e) {
            file >> g >> etype;
            std::getline(file, line);  // finish line
            ++nb_etype[etype];
            elements_max_glb_idx = std::max(elements_max_glb_idx, g);
        }

        // Allocate data structures for quads, triags, edges

        int nb_quads  = nb_etype[QUAD];
        int nb_triags = nb_etype[TRIAG];
        //int nb_edges  = nb_etype[LINE];

        mesh::Elements& quads  = mesh.cells().elements(mesh.cells().add(make_element_type(QUAD), nb_quads));
        mesh::Elements& triags = mesh.cells().elements(mesh.cells().add(make_element_type(TRIAG), nb_triags));
        //        mesh::Elements& edges  = mesh.edges().elements( mesh.edges().add( make_element_type( LINE ), nb_edges ) );

        mesh::BlockConnectivity& quad_nodes  = quads.node_connectivity();
        mesh::BlockConnectivity& triag_nodes = triags.node_connectivity();
        //        mesh::BlockConnectivity& edge_nodes  = edges.node_connectivity();

        array::ArrayView<gidx_t, 1> quad_glb_idx = array::make_view<gidx_t, 1>(quads.global_index());
        array::ArrayView<int, 1> quad_part       = array::make_view<int, 1>(quads.partition());

        array::ArrayView<gidx_t, 1> triag_glb_idx = array::make_view<gidx_t, 1>(triags.global_index());
        array::ArrayView<int, 1> triag_part       = array::make_view<int, 1>(triags.partition());

        //        array::ArrayView<gidx_t, 1> edge_glb_idx = array::make_view<gidx_t, 1>( edges.global_index() );
        //        array::ArrayView<int, 1> edge_part       = array::make_view<int, 1>( edges.partition() );

        // Now read all elements
        file.seekg(position, std::ios::beg);
        int gn0, gn1, gn2, gn3;
        int quad = 0, triag = 0;
        int ntags, tags[100];
        for (int e = 0; e < nb_elements; ++e) {
            file >> g >> etype >> ntags;
            for (int t = 0; t < ntags; ++t) {
                file >> tags[t];
            }
            int part = 0;
            if (ntags > 3) {
                part = std::max(part, *std::max_element(tags + 3, tags + ntags - 1));  // one positive, others negative
            }

            idx_t enodes[4] = {-1, -1, -1, -1};

            switch (etype) {
                case (QUAD):
                    file >> gn0 >> gn1 >> gn2 >> gn3;
                    quad_glb_idx(quad) = g;
                    quad_part(quad)    = part;
                    enodes[0]          = glb_to_loc[gn0];
                    enodes[1]          = glb_to_loc[gn1];
                    enodes[2]          = glb_to_loc[gn2];
                    enodes[3]          = glb_to_loc[gn3];
                    quad_nodes.set(quad, enodes);
                    ++quad;
                    break;
                case (TRIAG):
                    file >> gn0 >> gn1 >> gn2;
                    triag_glb_idx(triag) = g;
                    triag_part(triag)    = part;
                    enodes[0]            = glb_to_loc[gn0];
                    enodes[1]            = glb_to_loc[gn1];
                    enodes[2]            = glb_to_loc[gn2];
                    triag_nodes.set(triag, enodes);
                    ++triag;
                    break;
                case (LINE):
                    file >> gn0 >> gn1;
                    //                    edge_glb_idx( edge ) = g;
                    //                    edge_part( edge )    = part;
                    //                    enodes[0]            = glb_to_loc[gn0];
                    //                    enodes[1]            = glb_to_loc[gn1];
                    //                    edge_nodes.set( edge, enodes );
                    //                    ++edge;
                    break;
                case (POINT):
                    file >> gn0;
                    break;
                default:
                    std::cout << "etype " << etype << std::endl;
                    throw_Exception("ERROR: element type not supported", Here());
            }
        }
    }

    file.close();
}

void GmshIO::write(const Mesh& mesh, const PathName& file_path) const {
    mpi::Scope scope(mesh.mpi_comm());
    int part = mesh.metadata().has("part") ? mesh.metadata().get<size_t>("part") : mpi::rank();

    std::string nodes_field  = options.get<std::string>("nodes");
    const mesh::Nodes& nodes = mesh.nodes();

    const Field coords_field = nodes.field(nodes_field);
    array::ArrayT<double> dummy_double(1, 1);
    array::ArrayT<idx_t> dummy_idx(1, 1);
    bool coords_is_idx    = coords_field.datatype().kind() == array::make_datatype<idx_t>().kind();
    auto coords           = array::make_view<const double, 2>(coords_is_idx ? dummy_double : coords_field.array());
    auto coords_idx       = array::make_view<const idx_t, 2>(coords_is_idx ? coords_field.array() : dummy_idx);
    bool coords_is_lonlat = coords_field.name() == "lonlat";

    auto glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());

    const idx_t surfdim = nodes.field(nodes_field).shape(1);  // nb of variables in coords

    bool include_patch = (surfdim == 3);
    NodeInclusionPolicy node_inclusion(nodes);


    ATLAS_ASSERT(surfdim == 2 || surfdim == 3);

    Log::debug() << "writing mesh to gmsh file " << file_path << std::endl;

    const bool binary                      = !options.get<bool>("ascii");
    const bool element_partition_as_entity = options.get<bool>("element_partition_as_entity");

    gidx_t local_max_cell_tag = 0;
    auto cell_global_index    = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
        local_max_cell_tag = std::max(local_max_cell_tag, cell_global_index(cell));
    }
    gidx_t max_cell_tag = 0;
    mpi::comm().allReduce(local_max_cell_tag, max_cell_tag, eckit::mpi::max());

    gidx_t edge_tag_offset = 10;
    while (edge_tag_offset <= max_cell_tag) {
        if (edge_tag_offset > std::numeric_limits<int>::max() / 10) {
            ATLAS_THROW_EXCEPTION("No 32-bit Gmsh element tag range remains above the cell tags");
        }
        edge_tag_offset *= 10;
    }
    Field edge_element_tags = gmsh_edge_tags(mesh, edge_tag_offset);

    ATLAS_DEBUG_VAR(element_partition_as_entity);
    ATLAS_DEBUG_VAR(binary);

    openmode mode = std::ios::out;
    if (binary) {
        mode = std::ios::out | std::ios::binary;
    }
    GmshFile file(file_path, mode, part);

    // Header
    if (binary) {
        write_header_binary(file);
    }
    else {
        write_header_ascii(file);
    }

    struct ElementBlock {
        const mesh::Elements* elements;
        const Field* global_index;
        int dimension;
        int gmsh_type;
        size_t nb_nodes;
        std::vector<idx_t> included;
        std::map<int, std::vector<idx_t>> ghosts_by_owner;
        std::map<int, std::vector<idx_t>> by_owner;
    };

    std::vector<ElementBlock> element_blocks;
    std::vector<const mesh::HybridElements*> grouped_elements;
    if (options.get<bool>("elements")) {
        grouped_elements.push_back(&mesh.cells());
    }
    if (options.get<bool>("edges")) {
        grouped_elements.push_back(&mesh.edges());
    }

    for (const mesh::HybridElements* hybrid : grouped_elements) {
        for (idx_t etype = 0; etype < hybrid->nb_types(); ++etype) {
            const mesh::Elements& elements        = hybrid->elements(etype);
            const mesh::ElementType& element_type = elements.element_type();
            ElementInclusionPolicy element_inclusion(options, part, include_patch, coords_is_lonlat, coords, elements);
            const Field& element_global_index = hybrid == &mesh.edges() ? edge_element_tags : hybrid->global_index();
            ElementBlock block{&elements, &element_global_index, 2, 0, element_inclusion.nb_nodes(), {}, {}, {}};
            auto element_owner = elements.view<int, 1>(elements.partition());
            if (element_type.name() == "Line") {
                block.dimension = 1;
                block.gmsh_type = LINE;
            }
            else if (element_type.name() == "Triangle") {
                block.gmsh_type = TRIAG;
            }
            else if (element_type.name() == "Quadrilateral") {
                block.gmsh_type = QUAD;
            }
            else if (element_type.name() == "Pentagon") {
                block.gmsh_type = QUAD;
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }

            std::vector<ElementInclusion> inclusion(elements.size());
            std::map<int, size_t> owner_counts;
            size_t included_count = 0;
            for (idx_t elem = 0; elem < elements.size(); ++elem) {
                inclusion[elem] = element_inclusion(elem);
                switch (inclusion[elem]) {
                    case ElementInclusion::OWNED:
                        element_partition_as_entity ? ++owner_counts[element_owner(elem)] : ++included_count;
                        break;
                    case ElementInclusion::GHOST:
                        ++owner_counts[element_owner(elem)];
                        break;
                    case ElementInclusion::SKIP:
                        break;
                }
            }

            if (element_partition_as_entity) {
                for (const auto& [owner, count] : owner_counts) {
                    block.by_owner[owner].reserve(count);
                }
            }
            else {
                block.included.reserve(included_count);
                for (const auto& [owner, count] : owner_counts) {
                    block.ghosts_by_owner[owner].reserve(count);
                }
            }
            for (idx_t elem = 0; elem < elements.size(); ++elem) {
                switch (inclusion[elem]) {
                    case ElementInclusion::OWNED:
                        (element_partition_as_entity ? block.by_owner[element_owner(elem)] : block.included)
                            .push_back(elem);
                        break;
                    case ElementInclusion::GHOST:
                        (element_partition_as_entity ? block.by_owner[element_owner(elem)]
                                                     : block.ghosts_by_owner[element_owner(elem)])
                            .push_back(elem);
                        break;
                    case ElementInclusion::SKIP:
                        break;
                }
            }
            if (element_partition_as_entity ? !block.by_owner.empty()
                                            : !block.included.empty() || !block.ghosts_by_owner.empty()) {
                element_blocks.emplace_back(std::move(block));
            }
        }
    }

    const idx_t nb_nodes                     = static_cast<idx_t>(node_inclusion.size());
    const int partition_tag                  = part;
    const int canonical_entity_tag           = part + 2;
    const int node_entity_tag = element_partition_as_entity ? partition_tag : canonical_entity_tag;
    const int nb_partitions =
        std::max(part + 1, mesh.metadata().has("nb_parts") ? mesh.metadata().get<int>("nb_parts")
                                                           : static_cast<int>(mpi::size()));
    const bool has_curve_blocks = std::any_of(element_blocks.begin(), element_blocks.end(),
                                              [](const auto& block) { return block.dimension == 1; });
    std::map<int, bool> curve_entity_owners;
    std::map<int, bool> surface_entity_owners;
    if (element_partition_as_entity) {
        if (nb_nodes) {
            surface_entity_owners[partition_tag] = true;
        }
        for (const auto& block : element_blocks) {
            auto& owners = block.dimension == 1 ? curve_entity_owners : surface_entity_owners;
            for ([[maybe_unused]] const auto& [owner, elements] : block.by_owner) {
                owners[owner] = true;
            }
        }
    }
    else {
        if (nb_nodes) {
            surface_entity_owners[part] = true;
        }
        for (const auto& block : element_blocks) {
            auto& owners = block.dimension == 1 ? curve_entity_owners : surface_entity_owners;
            if (!block.included.empty()) {
                owners[part] = true;
            }
            for ([[maybe_unused]] const auto& [owner, elements] : block.ghosts_by_owner) {
                owners[owner] = true;
            }
        }
    }
    const size_t nb_ghost_entities = 0;

    double xyz_min[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max()};
    double xyz_max[3] = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                         -std::numeric_limits<double>::max()};
    auto node_xyz     = [&](idx_t node) {
        std::array<double, 3> xyz{0., 0., 0.};
        for (idx_t dimension = 0; dimension < surfdim; ++dimension) {
            xyz[dimension] = coords_is_idx ? coords_idx(node, dimension) : coords(node, dimension);
        }
        return xyz;
    };
    for (idx_t node : node_inclusion.included()) {
        auto xyz = node_xyz(node);
        for (idx_t dimension = 0; dimension < 3; ++dimension) {
            xyz_min[dimension] = std::min(xyz_min[dimension], xyz[dimension]);
            xyz_max[dimension] = std::max(xyz_max[dimension], xyz[dimension]);
        }
    }
    if (nb_nodes == 0) {
        std::fill(xyz_min, xyz_min + 3, 0.);
        std::fill(xyz_max, xyz_max + 3, 0.);
    }

    auto write_bounding_box = [&]() {
        file << xyz_min[0] << " " << xyz_min[1] << " " << xyz_min[2] << " " << xyz_max[0] << " " << xyz_max[1] << " "
             << xyz_max[2];
    };
    auto write_binary = [&](const auto& value) { file.write(reinterpret_cast<const char*>(&value), sizeof(value)); };
    auto write_binary_bounding_box = [&]() {
        file.write(reinterpret_cast<const char*>(xyz_min), sizeof(xyz_min));
        file.write(reinterpret_cast<const char*>(xyz_max), sizeof(xyz_max));
    };

    file << "$Entities\n";
    if (binary) {
        const size_t entity_counts[4] = {
            0, element_partition_as_entity ? curve_entity_owners.size() : (has_curve_blocks ? 1ul : 0ul),
            element_partition_as_entity ? surface_entity_owners.size() : 1ul, 0};
        file.write(reinterpret_cast<const char*>(entity_counts), sizeof(entity_counts));
        auto write_entity_record = [&](int tag) {
            const size_t nb_physical   = 1;
            const int physical_tag     = 1;
            const size_t nb_boundaries = 0;
            write_binary(tag);
            write_binary_bounding_box();
            write_binary(nb_physical);
            write_binary(physical_tag);
            write_binary(nb_boundaries);
        };
        if (element_partition_as_entity) {
            for ([[maybe_unused]] const auto& [owner, present] : curve_entity_owners) {
                write_entity_record(owner);
            }
            for ([[maybe_unused]] const auto& [owner, present] : surface_entity_owners) {
                write_entity_record(owner);
            }
        }
        else {
            if (has_curve_blocks) {
                write_entity_record(1);  // parent curve entity
            }
            write_entity_record(1);  // parent surface entity
        }
        file << "\n";
    }
    else {
        file << "0 " << (element_partition_as_entity ? curve_entity_owners.size() : (has_curve_blocks ? 1 : 0)) << " "
             << (element_partition_as_entity ? surface_entity_owners.size() : 1) << " 0\n";
        auto write_entity_record = [&](int tag) {
            file << tag << " ";
            write_bounding_box();
            file << " 1 1 0\n";
        };
        if (element_partition_as_entity) {
            for ([[maybe_unused]] const auto& [owner, present] : curve_entity_owners) {
                write_entity_record(owner);
            }
            for ([[maybe_unused]] const auto& [owner, present] : surface_entity_owners) {
                write_entity_record(owner);
            }
        }
        else {
            if (has_curve_blocks) {
                write_entity_record(1);
            }
            write_entity_record(1);
        }
    }
    file << "$EndEntities\n";

    constexpr int dim_1d = 1;
    constexpr int dim_2d = 2;
    if (!element_partition_as_entity) {
        file << "$PartitionedEntities\n";
        if (binary) {
            const size_t partitions       = nb_partitions;
            const size_t entity_counts[4] = {0, curve_entity_owners.size(), surface_entity_owners.size(), 0};
            write_binary(partitions);
            write_binary(nb_ghost_entities);
            file.write(reinterpret_cast<const char*>(entity_counts), sizeof(entity_counts));
            auto write_partitioned_entity = [&](int owner, int dimension) {
                const int parent_tag_value        = 1;
                const size_t entity_nb_partitions = 1;
                const size_t nb_physical          = 1;
                const int physical_tag            = 1;
                const size_t nb_boundaries        = 0;
                write_binary(owner + 2);
                write_binary(dimension);
                write_binary(parent_tag_value);
                write_binary(entity_nb_partitions);
                write_binary(owner);
                write_binary_bounding_box();
                write_binary(nb_physical);
                write_binary(physical_tag);
                write_binary(nb_boundaries);
            };
            for ([[maybe_unused]] const auto& [owner, present] : curve_entity_owners) {
                write_partitioned_entity(owner, dim_1d);
            }
            for ([[maybe_unused]] const auto& [owner, present] : surface_entity_owners) {
                write_partitioned_entity(owner, dim_2d);
            }
            file << "\n";
        }
        else {
            file << nb_partitions << "\n";
            file << nb_ghost_entities << "\n";
            file << "0 " << curve_entity_owners.size() << " " << surface_entity_owners.size() << " 0\n";
            auto write_partitioned_entity = [&](int owner, int dimension) {
                file << owner + 2 << " " << dimension << " 1 1 " << owner << " ";
                write_bounding_box();
                file << " 1 1 0\n";
            };
            for ([[maybe_unused]] const auto& [owner, present] : curve_entity_owners) {
                write_partitioned_entity(owner, dim_1d);
            }
            for ([[maybe_unused]] const auto& [owner, present] : surface_entity_owners) {
                write_partitioned_entity(owner, dim_2d);
            }
        }
        file << "$EndPartitionedEntities\n";
    }

    gidx_t min_node_tag = node_inclusion.empty() ? 0 : glb_idx(node_inclusion.included().front());
    gidx_t max_node_tag = min_node_tag;
    for (idx_t node : node_inclusion.included()) {
        min_node_tag = std::min(min_node_tag, glb_idx(node));
        max_node_tag = std::max(max_node_tag, glb_idx(node));
    }
    file << "$Nodes\n";
    if (binary) {
        const size_t node_header[4] = {nb_nodes ? 1ul : 0ul, static_cast<size_t>(nb_nodes),
                                       static_cast<size_t>(min_node_tag), static_cast<size_t>(max_node_tag)};
        file.write(reinterpret_cast<const char*>(node_header), sizeof(node_header));
        if (nb_nodes) {
            const int dimension  = 2;
            const int parametric = 0;
            const size_t count   = nb_nodes;
            write_binary(dimension);
            write_binary(node_entity_tag);
            write_binary(parametric);
            write_binary(count);
            for (idx_t node : node_inclusion.included()) {
                const size_t tag = glb_idx(node);
                write_binary(tag);
            }
            for (idx_t node : node_inclusion.included()) {
                auto xyz = node_xyz(node);
                file.write(reinterpret_cast<const char*>(xyz.data()), sizeof(double) * xyz.size());
            }
        }
        file << "\n";
    }
    else {
        file << (nb_nodes ? 1 : 0) << " " << nb_nodes << " " << min_node_tag << " " << max_node_tag << "\n";
        if (nb_nodes) {
            file << "2 " << node_entity_tag << " 0 " << nb_nodes << "\n";
            for (idx_t node : node_inclusion.included()) {
                file << glb_idx(node) << "\n";
            }
            for (idx_t node : node_inclusion.included()) {
                auto xyz = node_xyz(node);
                file << xyz[0] << " " << xyz[1] << " " << xyz[2] << "\n";
            }
        }
    }
    file << "$EndNodes\n";

    size_t nb_elements     = 0;
    gidx_t min_element_tag = 0;
    gidx_t max_element_tag = 0;
    for (const auto& block : element_blocks) {
        auto elems_glb_idx = block.elements->view<gidx_t, 1>(*block.global_index);
        auto update_element_range = [&](idx_t elem) {
            const gidx_t tag = elems_glb_idx(elem);
            if (nb_elements++ == 0) {
                min_element_tag = tag;
                max_element_tag = tag;
            }
            else {
                min_element_tag = std::min(min_element_tag, tag);
                max_element_tag = std::max(max_element_tag, tag);
            }
        };
        if (element_partition_as_entity) {
            for (const auto& [owner, elements] : block.by_owner) {
                for (idx_t elem : elements) {
                    update_element_range(elem);
                }
            }
        }
        else {
            for (idx_t elem : block.included) {
                update_element_range(elem);
            }
            for (const auto& [owner, elements] : block.ghosts_by_owner) {
                for (idx_t elem : elements) {
                    update_element_range(elem);
                }
            }
        }
    }
    const size_t nb_element_blocks = std::accumulate(
        element_blocks.begin(), element_blocks.end(), size_t{0}, [&](size_t count, const auto& block) {
            return count + (element_partition_as_entity ? block.by_owner.size()
                                                        : !block.included.empty() + block.ghosts_by_owner.size());
        });
    file << "$Elements\n";
    if (binary) {
        const size_t element_header[4] = {nb_element_blocks, nb_elements, static_cast<size_t>(min_element_tag),
                                          static_cast<size_t>(max_element_tag)};
        file.write(reinterpret_cast<const char*>(element_header), sizeof(element_header));
        for (const auto& block : element_blocks) {
            const auto& node_connectivity = block.elements->node_connectivity();
            auto elems_glb_idx            = block.elements->view<gidx_t, 1>(*block.global_index);
            auto write_block              = [&](const std::vector<idx_t>& elements, int block_entity_tag) {
                if (elements.empty()) {
                    return;
                }
                const size_t count = elements.size();
                write_binary(block.dimension);
                write_binary(block_entity_tag);
                write_binary(block.gmsh_type);
                write_binary(count);
                for (idx_t elem : elements) {
                    const size_t element_tag = elems_glb_idx(elem);
                    write_binary(element_tag);
                    for (idx_t node = 0; node < block.nb_nodes; ++node) {
                        const size_t node_tag = glb_idx(node_connectivity(elem, node));
                        write_binary(node_tag);
                    }
                }
            };
            if (element_partition_as_entity) {
                for (const auto& [owner, elements] : block.by_owner) {
                    write_block(elements, owner);
                }
            }
            else {
                write_block(block.included, canonical_entity_tag);
                for (const auto& [owner, elements] : block.ghosts_by_owner) {
                    write_block(elements, owner + 2);
                }
            }
        }
        file << "\n";
    }
    else {
        file << nb_element_blocks << " " << nb_elements << " " << min_element_tag << " " << max_element_tag << "\n";
        for (const auto& block : element_blocks) {
            const auto& node_connectivity = block.elements->node_connectivity();
            auto elems_glb_idx            = block.elements->view<gidx_t, 1>(*block.global_index);
            auto write_block              = [&](const std::vector<idx_t>& elements, int block_entity_tag) {
                if (elements.empty()) {
                    return;
                }
                file << block.dimension << " " << block_entity_tag << " " << block.gmsh_type << " "
                     << elements.size() << "\n";
                for (idx_t elem : elements) {
                    file << elems_glb_idx(elem);
                    for (idx_t node = 0; node < block.nb_nodes; ++node) {
                        file << " " << glb_idx(node_connectivity(elem, node));
                    }
                    file << "\n";
                }
            };
            if (element_partition_as_entity) {
                for (const auto& [owner, elements] : block.by_owner) {
                    write_block(elements, owner);
                }
            }
            else {
                write_block(block.included, canonical_entity_tag);
                for (const auto& [owner, elements] : block.ghosts_by_owner) {
                    write_block(elements, owner + 2);
                }
            }
        }
    }
    file << "$EndElements\n";

    size_t nb_ghost_elements = 0;
    for (const auto& block : element_blocks) {
        for (const auto& [owner, elements] : block.ghosts_by_owner) {
            nb_ghost_elements += elements.size();
        }
    }
    if (nb_ghost_elements && !element_partition_as_entity) {
        file << "$GhostElements\n";
        if (binary) {
            write_binary(nb_ghost_elements);
            for (const auto& block : element_blocks) {
                auto elems_glb_idx = block.elements->view<gidx_t, 1>(*block.global_index);
                for (const auto& [owner, elements] : block.ghosts_by_owner) {
                    for (idx_t elem : elements) {
                        const size_t element_tag         = elems_glb_idx(elem);
                        const size_t nb_ghost_partitions = 1;
                        write_binary(element_tag);
                        write_binary(owner);
                        write_binary(nb_ghost_partitions);
                        write_binary(partition_tag);
                    }
                }
            }
            file << "\n";
        }
        else {
            file << nb_ghost_elements << "\n";
            for (const auto& block : element_blocks) {
                auto elems_glb_idx = block.elements->view<gidx_t, 1>(*block.global_index);
                for (const auto& [owner, elements] : block.ghosts_by_owner) {
                    for (idx_t elem : elements) {
                        file << elems_glb_idx(elem) << " " << owner << " 1 " << partition_tag << "\n";
                    }
                }
            }
        }
        file << "$EndGhostElements\n";
    }
    file << std::flush;

    // Optional mesh information file
    if (options.has("info") && options.get<bool>("info")) {
        PathName mesh_info(file_path);
        mesh_info = mesh_info.dirName() + "/" + mesh_info.baseName(false) + "_info.msh";

        const std::vector<std::string> extra_fields = {"partition", "water", "dual_volumes", "dual_delta_sph",
                                                       "ghost",     "halo",  "remote_index"};
        {
            functionspace::NodeColumns function_space(mesh);
            FieldNamePrefixScope field_name_prefix_scope("nodes.");
            FieldSet fieldset;
            auto lat =
                array::make_view<double, 1>(fieldset.add(Field("lat", array::make_datatype<double>(), {nodes.size()})));
            auto lon =
                array::make_view<double, 1>(fieldset.add(Field("lon", array::make_datatype<double>(), {nodes.size()})));

            auto lonlat = array::make_view<double, 2>(nodes.lonlat());
            for (idx_t n = 0; n < nodes.size(); ++n) {
                lon(n) = lonlat(n, 0);
                lat(n) = lonlat(n, 1);
            }
            write(fieldset, function_space, mesh_info, std::ios_base::out);
            for (const auto& field_name : extra_fields) {
                if (nodes.has_field(field_name)) {
                    write(nodes.field(field_name), function_space, mesh_info, std::ios_base::app);
                }
            }
        }

        {
            functionspace::CellColumns function_space(mesh);
            FieldNamePrefixScope field_name_prefix_scope("cells.");
            for (const auto& field_name : extra_fields) {
                if (mesh.cells().has_field(field_name)) {
                    write(mesh.cells().field(field_name), function_space, mesh_info, std::ios_base::app);
                }
            }
        }

        if (mesh.edges().size()) {
            functionspace::EdgeColumns function_space(mesh);
            FieldNamePrefixScope field_name_prefix_scope("edges.");
            ElementGlobalIndexScope element_global_index_scope(edge_element_tags);
            for (const auto& field_name : extra_fields) {
                if (mesh.edges().has_field(field_name)) {
                    write(mesh.edges().field(field_name), function_space, mesh_info, std::ios_base::app);
                }
            }
        }
    }
    file.close();
}

// ----------------------------------------------------------------------------
void GmshIO::write(const Field& field, const PathName& file_path, openmode mode) const {
    if (!field.functionspace()) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }

    else if (functionspace::NodeColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else if (functionspace::StructuredColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else if (functionspace::CellColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else if (functionspace::EdgeColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else {
        std::stringstream msg;
        msg << "Field [" << field.name() << "] has functionspace [" << field.functionspace().type()
            << "] but requires a [functionspace::NodeColumns "
            << "or functionspace::StructuredColumns or functionspace::CellColumns "
            << "or functionspace::EdgeColumns]";

        throw_AssertionFailed(msg.str(), Here());
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const Field& field, const functionspace::NodeColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------

void GmshIO::write_delegate(const Field& field, const functionspace::NoFunctionSpace& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------

void GmshIO::write_delegate(const Field& field, const functionspace::CellColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const Field& field, const functionspace::EdgeColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const Field& field, const functionspace::StructuredColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::NodeColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        binary ? write_header_binary(file) : write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}

void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::NoFunctionSpace& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        binary ? write_header_binary(file) : write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}

void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::CellColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        binary ? write_header_binary(file) : write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_elems<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_elems<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_elems<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_elems<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}
// ----------------------------------------------------------------------------

void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::EdgeColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    if (is_new_file) {
        binary ? write_header_binary(file) : write_header_ascii(file);
    }

    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_elems<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_elems<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_elems<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_elems<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::StructuredColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));

    if (binary) {
        mode |= std::ios_base::binary;
    }

    bool gather = options.has("gather") ? options.get<bool>("gather") : false;

    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        binary ? write_header_binary(file) : write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }

    file.close();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write(const FieldSet& fieldset, const FunctionSpace& funcspace, const eckit::PathName& file_path,
                   openmode mode) const {
    mpi::Scope scope(funcspace.mpi_comm());
    if (functionspace::NodeColumns(funcspace)) {
        write_delegate(fieldset, functionspace::NodeColumns(funcspace), file_path, mode);
    }
    else if (functionspace::StructuredColumns(funcspace)) {
        write_delegate(fieldset, functionspace::StructuredColumns(funcspace), file_path, mode);
    }
    else if (functionspace::CellColumns(funcspace)) {
        write_delegate(fieldset, functionspace::CellColumns(funcspace), file_path, mode);
    }
    else if (functionspace::EdgeColumns(funcspace)) {
        write_delegate(fieldset, functionspace::EdgeColumns(funcspace), file_path, mode);
    }
    else if (not funcspace) {
        write_delegate(fieldset, functionspace::NoFunctionSpace(), file_path, mode);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write(const Field& field, const FunctionSpace& funcspace, const eckit::PathName& file_path,
                   openmode mode) const {
    mpi::Scope scope(funcspace.mpi_comm());
    if (functionspace::NodeColumns(funcspace)) {
        write_delegate(field, functionspace::NodeColumns(funcspace), file_path, mode);
    }
    else if (functionspace::StructuredColumns(funcspace)) {
        write_delegate(field, functionspace::StructuredColumns(funcspace), file_path, mode);
    }
    else if (functionspace::CellColumns(funcspace)) {
        write_delegate(field, functionspace::CellColumns(funcspace), file_path, mode);
    }
    else if (functionspace::EdgeColumns(funcspace)) {
        write_delegate(field, functionspace::EdgeColumns(funcspace), file_path, mode);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}
// ----------------------------------------------------------------------------

class GmshFortranInterface {
public:
    static Mesh::Implementation* atlas__Gmsh__read(GmshIO* This, char* file_path);
    static void atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path);
    static Mesh::Implementation* atlas__read_gmsh(char* file_path);
    static void atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path);
    static void atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset,
                                           functionspace::FunctionSpaceImpl* functionspace, char* file_path, int mode);
    static void atlas__write_gmsh_field(const field::FieldImpl* field, functionspace::FunctionSpaceImpl* functionspace,
                                        char* file_path, int mode);
};

Mesh::Implementation* GmshFortranInterface::atlas__Gmsh__read(GmshIO* This, char* file_path) {
    Mesh::Implementation* m;
    {
        Mesh mesh = This->read(PathName(file_path));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

void GmshFortranInterface::atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path) {
    Mesh m(mesh);
    This->write(m, PathName(file_path));
}

Mesh::Implementation* GmshFortranInterface::atlas__read_gmsh(char* file_path) {
    Mesh::Implementation* m;
    {
        Mesh mesh = GmshIO().read(PathName(file_path));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

void GmshFortranInterface::atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path) {
    GmshIO writer;
    writer.write(mesh, PathName(file_path));
}

void GmshFortranInterface::atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset,
                                                      functionspace::FunctionSpaceImpl* functionspace, char* file_path,
                                                      int /*mode*/) {
    GmshIO writer;
    writer.write(fieldset, functionspace, PathName(file_path));
}

void GmshFortranInterface::atlas__write_gmsh_field(const field::FieldImpl* field,
                                                   functionspace::FunctionSpaceImpl* functionspace, char* file_path,
                                                   int /*mode*/) {
    GmshIO writer;
    writer.write(field, functionspace, PathName(file_path));
}

extern "C" {

// ----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
// ----------------------------------------------------------------------------
GmshIO* atlas__Gmsh__new() {
    return new GmshIO();
}

void atlas__Gmsh__delete(GmshIO* This) {
    delete This;
}

Mesh::Implementation* atlas__Gmsh__read(GmshIO* This, char* file_path) {
    return GmshFortranInterface::atlas__Gmsh__read(This, file_path);
}

void atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path) {
    GmshFortranInterface::atlas__Gmsh__write(This, mesh, file_path);
}

Mesh::Implementation* atlas__read_gmsh(char* file_path) {
    return GmshFortranInterface::atlas__read_gmsh(file_path);
}

void atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path) {
    GmshFortranInterface::atlas__write_gmsh_mesh(mesh, file_path);
}

void atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset, functionspace::FunctionSpaceImpl* functionspace,
                                char* file_path, int mode) {
    GmshFortranInterface::atlas__write_gmsh_fieldset(fieldset, functionspace, file_path, mode);
}

void atlas__write_gmsh_field(const field::FieldImpl* field, functionspace::FunctionSpaceImpl* functionspace,
                             char* file_path, int mode) {
    GmshFortranInterface::atlas__write_gmsh_field(field, functionspace, file_path, mode);
}
}
// ----------------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
