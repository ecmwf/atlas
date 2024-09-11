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
#include "atlas/field/MultiFieldCreatorIFS.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

MultiFieldCreatorIFS::MultiFieldCreatorIFS() {}

MultiFieldCreatorIFS::MultiFieldCreatorIFS(const eckit::Configuration& config) {}

MultiFieldCreatorIFS::~MultiFieldCreatorIFS() = default;

MultiFieldImpl* MultiFieldCreatorIFS::create(const array::DataType datatype, const std::vector<int>& shape,
        const std::vector<std::string>& var_names) const {
    ATLAS_NOTIMPLEMENTED;
    return nullptr;
}

/*
template<typename T>
MultiFieldImpl* MultiFieldCreatorIFS::create(const std::vector<int> shape, const std::vector<std::string> var_names) const {
    const int dim = shape.size();
    const int nvar = var_names.size();
    ATLAS_ASSERT(nvar > 0);

    ATLAS_ASSERT(nvar > 0 && dim > 2, "MultiField must have at least one field name.");

    int varidx = -1;
    for (int i = 0; i < dim; i++) {
        if (shape[i] == -1) {
            varidx = i;
            break;
        }
    }
    int nlev = 0;
    if (varidx > 0) {
        nlev = shape[2];
    }

    const int nblk = shape[0];
    const int nproma = shape[dim - 1];

    util::Config config;
    config.set("type", "MultiFieldCreatorIFS");
    config.set("datatype", array::make_datatype<T>().str());
    config.set("nproma", nproma);
    config.set("nblk", nblk);
    config.set("nlev", nlev);
    config.set("ngptot", nblk * nproma);

    std::vector<util::Config> fconfigs(nvar);
    for (int i = 0; i < nvar; i++) {
        fconfigs[i].set("name", var_names[i]);
    }
    config.set("fields", fconfigs);
    return create(config);
}
*/

MultiFieldImpl* MultiFieldCreatorIFS::create(const eckit::Configuration& config) const {
    long nproma;
    config.get("nproma", nproma);
    long nlev;
    config.get("nlev", nlev);
    long nblk = 0;
    if (config.has("nblk")) {
        config.get("nblk", nblk);
    }
    else if (config.has("ngptot")) {
        long ngptot;
        config.get("ngptot", ngptot);
        nblk        = std::ceil(static_cast<double>(ngptot) / static_cast<double>(nproma));
    }
    else {
        ATLAS_THROW_EXCEPTION("Configuration not found: ngptot or nblk");
    }
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

    auto fields = config.getSubConfigurations("fields");
    long nfld   = 0;
    for (const auto& field_params : fields) {
        long nvar = 1;
        field_params.get("nvar", nvar);
        nfld += nvar;
    }
    array::ArrayShape multiarray_shape = ( (nlev > 0 and nfld > 0) ? array::make_shape(nblk, nfld, nlev, nproma) : 
        ( (nlev > 0) ? array::make_shape(nblk, nlev, nproma) : ( (nfld > 0) ? array::make_shape(nblk, nfld, nproma) : 
                                                                 array::make_shape(nblk, nproma) ) ) );

    MultiFieldImpl* multifield = new MultiFieldImpl{array::ArraySpec{datatype, multiarray_shape}, config};
    auto& multiarray = multifield->array();

    size_t multiarray_field_idx = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        Field field;
        size_t field_vars = 1;
        if (fields[i].get("nvar", field_vars)) {
            array::ArrayShape field_shape =
                ( (nlev > 0 and field_vars > 0) ?
                  array::make_shape(multiarray.shape(0), field_vars, multiarray.shape(2), multiarray.shape(3)) :
                  ( nlev > 0 ?
                  array::make_shape(multiarray.shape(0), multiarray.shape(1), multiarray.shape(2)) :
                  ( field_vars > 0 ?
                  array::make_shape(multiarray.shape(0), field_vars, multiarray.shape(2)) :
                  array::make_shape(multiarray.shape(0), multiarray.shape(1)) ) ) );
            array::ArrayShape multiarray_shape =
                ( (nlev > 0 and field_vars > 0) ? array::make_shape(nblk, field_vars, nlev, nproma) : 
                ( (nlev > 0) ? array::make_shape(nblk, nlev, nproma) : ( (field_vars > 0) ?
                array::make_shape(nblk, field_vars, nproma) : array::make_shape(nblk, nproma) ) ) );
            auto field_strides    = multiarray.strides();
            auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            constexpr auto all = array::Range::all();
            const auto range   = array::Range(multiarray_field_idx, multiarray_field_idx + field_vars);
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                if (nlev > 0 and field_vars > 0) {
                    auto slice = array::make_view<double, 4>(multiarray).slice(all, range, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (nlev > 0) {
                    auto slice = array::make_view<double, 3>(multiarray).slice(all, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (field_vars > 0) {
                    auto slice = array::make_view<double, 3>(multiarray).slice(all, range, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else {
                    auto slice = array::make_view<double, 2>(multiarray).slice(all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                if (nlev > 0 and field_vars > 0) {
                    auto slice = array::make_view<float, 4>(multiarray).slice(all, range, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (nlev > 0) {
                    auto slice = array::make_view<float, 3>(multiarray).slice(all, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (field_vars > 0) {
                    auto slice = array::make_view<float, 3>(multiarray).slice(all, range, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else {
                    auto slice = array::make_view<float, 2>(multiarray).slice(all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            field.set_variables(field_vars);
        }
        else {
            array::ArraySpec field_array_spec;
            if (nlev > 0) {
                auto field_shape   = array::make_shape(multiarray.shape(0), multiarray.shape(2), multiarray.shape(3));
                auto field_strides = array::make_strides(multiarray.stride(0), multiarray.stride(2), multiarray.stride(3));
                field_array_spec = array::ArraySpec(field_shape, field_strides);
            }
            else if (field_vars > 0) {
                auto field_shape = array::make_shape(multiarray.shape(0), multiarray.shape(2));
                auto field_strides = array::make_strides(multiarray.stride(0), multiarray.stride(2));
                field_array_spec = array::ArraySpec(field_shape, field_strides);
            }

            constexpr auto all = array::Range::all();
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                if (nlev > 0 and field_vars > 0) {
                    auto slice = array::make_view<double, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (nlev > 0) {
                    auto slice = array::make_view<double, 3>(multiarray).slice(all, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (field_vars > 0) {
                    auto slice = array::make_view<double, 3>(multiarray).slice(all, multiarray_field_idx, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else {
                    auto slice = array::make_view<double, 2>(multiarray).slice(all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                if (nlev > 0 and field_vars > 0) {
                    auto slice = array::make_view<float, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (nlev > 0) {
                    auto slice = array::make_view<float, 3>(multiarray).slice(all, all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else if (field_vars > 0) {
                    auto slice = array::make_view<float, 3>(multiarray).slice(all, multiarray_field_idx, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
                else {
                    auto slice = array::make_view<float, 2>(multiarray).slice(all, all);
                    field      = Field(name, slice.data(), field_array_spec, config);
                }
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        field.set_levels(nlev);
        //field.set_blocks(nblk);

        multifield->add(field);

        multiarray_field_idx += field_vars;
    }
    std::string name;
    config.get("name", name);
    Log::debug() << "Creating IFS " << datatype.str() << " multifield: " << name << "[nblk=" << nblk << "][nvar=" << nfld
                 << "][nlev=" << nlev << "][nproma=" << nproma << "]\n";
    return multifield;
}

// ------------------------------------------------------------------

namespace {  
static MultiFieldCreatorBuilder<MultiFieldCreatorIFS> __MultiFieldCreatorIFS("MultiFieldCreatorIFS");
}

}  // namespace field
}  // namespace atlas
