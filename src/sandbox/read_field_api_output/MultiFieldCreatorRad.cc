/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Slavko Brdar
/// @date June 2024

#include <cmath>
#include <sstream>

#include "eckit/config/Parametrisation.h"

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/DataType.h"
#include "MultiFieldCreatorRad.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

static auto to_upper = [](const std::string& str) {
    std::string upper{str};
    for (auto & c: upper) c = std::toupper(c);
    return upper;
};


MultiFieldCreatorRad::MultiFieldCreatorRad() = default;

MultiFieldCreatorRad::MultiFieldCreatorRad(const eckit::Configuration& config) {}

MultiFieldCreatorRad::~MultiFieldCreatorRad() = default;

MultiFieldImpl* MultiFieldCreatorRad::create(const array::DataType datatype, const array::ArrayShape& shape,
        const std::vector<std::string>& var_names) const {
    ATLAS_NOTIMPLEMENTED;
    return nullptr;
}

MultiFieldImpl* MultiFieldCreatorRad::create(const eckit::Configuration& config) const {
    long dimension = config.getLong("dimension");
    long nproma = config.getLong("nproma");
    long nlev = config.getLong("nlev");
    long ngptot = config.getLong("ngptot");
    long nblk = std::ceil(static_cast<double>(ngptot) / static_cast<double>(nproma));
    long nlwemiss = config.getLong("nlwemiss", 0);
    long nsw = config.getLong("nsw", 0);
    long nrftotal_radgrid = config.getLong("nrftotal_radgrid",0);
    long nradaer = config.getLong("nradaer", 0);
    long nlwout = config.getLong("nlwout", 0);
    long nprogaer = config.getLong("nprogaer", 0);
    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (config.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        config.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    auto fields = config.getSubConfigurations("variables");

    auto is_number = [](const std::string & s) -> bool {
        for (char c : s) {
            if (!std::isdigit(c)) {
                return false;
            }
        }
        return true;
    };
    auto lookup_dim = [&](const std::string& dim) -> long {
        static std::map<std::string, long> map {
            {"nlev", nlev},
            {"nlev+1", nlev+1},
            {"nlwemiss", nlwemiss},
            {"nsw", nsw},
            {"nrftotal_radgrid",nrftotal_radgrid},
            {"nprogaer", nprogaer},
            {"nradaer",nradaer},
            {"nlwout",nlwout},
        };
        if (map.find(dim) == map.end()) {
            ATLAS_THROW_EXCEPTION("Could not match dimension '"<<dim<<"'");
        }
        return map.at(dim);
    };

    auto get_nvar = [&](const auto& field_params) {
        if (field_params.has("size")) {
            auto size = field_params.getLong("size");
            ATLAS_DEBUG_VAR(size);
            if (size == 0) {
                return size;
            }
            // return (field_params.getLong("size") / (nblk * nproma));
        }
        long nvar = 1;
        std::vector<std::string> dim_str;
        field_params.get("dim", dim_str);
        ATLAS_DEBUG_VAR(dim_str);
        for( const std::string& dim : dim_str ) {
            if (is_number(dim)) {
                nvar *= std::stol(dim);
            }
            else {
                nvar *= lookup_dim(dim);
            }
        }
        return nvar;
    };

    auto get_nfld = [&]() {
        long nfld   = 0;
        for (const auto& field : fields) {
            nfld += get_nvar(field);
        }
        return nfld;
    };
    long nfld = get_nfld();
    array::ArrayShape multiarray_shape = dimension == 3 ? array::make_shape(nblk, nfld, nproma)
                                       : array::make_shape(nblk, nproma);

    MultiFieldImpl* multifield = new MultiFieldImpl{array::ArraySpec{datatype, multiarray_shape}};
    auto& multiarray = multifield->array();

    size_t multiarray_field_idx = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        name = to_upper(name);
        Field field;
        ATLAS_DEBUG_VAR(name);
        long field_vars = get_nvar(fields[i]);
        ATLAS_DEBUG_VAR(field_vars);
        constexpr auto all = array::Range::all();
        if (field_vars == 0) {
            // empty field
            field      = Field(name, datatype, array::make_shape(multiarray.shape(0), 0, multiarray.shape(2)));
            field.set_levels(9999);
        }
        if (field_vars == 1) {
            // 2d field
            const auto field_shape      = array::make_shape(multiarray.shape(0), multiarray.shape(2));
            const auto field_strides    = array::make_strides(multiarray.stride(0), multiarray.stride(2));
            const auto field_array_spec = array::ArraySpec(field_shape, field_strides);
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 3>(multiarray).slice(all, multiarray_field_idx, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 3>(multiarray).slice(all, multiarray_field_idx, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        else {
            // 3d field
            const auto var_range     = array::Range(multiarray_field_idx, multiarray_field_idx + field_vars);
            const auto field_shape   = array::make_shape(multiarray.shape(0), field_vars, multiarray.shape(2));
            const auto field_strides = array::make_strides(multiarray.stride(0), multiarray.stride(1), multiarray.stride(2));
            const auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 3>(multiarray).slice(all, var_range, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 3>(multiarray).slice(all, var_range, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            field.set_levels(field_vars);
        }
        multifield->add(field);
        multiarray_field_idx += field_vars;
    }
    std::string name;
    config.get("name", name);
    Log::debug() << "Creating Radiation " << datatype.str() << " multifield: " << name << "[nblk=" << nblk << "][nfld=" << nfld
                 << "][nproma=" << nproma << "]\n";
    return multifield;
}

// ------------------------------------------------------------------

namespace {  
static MultiFieldCreatorBuilder<MultiFieldCreatorRad> __MultiFieldCreatorRad("MultiFieldCreatorRad");
}

}  // namespace field
}  // namespace atlas
