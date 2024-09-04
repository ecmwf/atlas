/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/FieldCreatorIFS.h"

#include <cmath>
#include <sstream>

#include "eckit/config/Parametrisation.h"

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/DataType.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

FieldImpl* FieldCreatorIFS::createField(const eckit::Parametrisation& params) const {
    size_t ngptot;
    size_t nblk;
    size_t nvar   = 1;
    size_t nproma = 1;
    size_t nlev   = 1;

    if (!params.get("ngptot", ngptot)) {
        throw_Exception("Could not find parameter 'ngptot' in Parametrisation");
    }
    params.get("nproma", nproma);
    params.get("nlev", nlev);
    params.get("nvar", nvar);

    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (params.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        params.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    nblk = std::ceil(static_cast<double>(ngptot + nproma - 1) / static_cast<double>(nproma));

    array::ArrayShape s;
    bool fortran(false);
    params.get("fortran", fortran);
    printf(" FieldCreatorIFS nlev, nvar: %zu, %zu", nlev, nvar);

    if (nlev > 0 && nvar > 0) {
        s = ( fortran ? array::make_shape(nproma, nlev, nvar, nblk) : array::make_shape(nblk, nvar, nlev, nproma) );
    }
    else if (nlev > 0) {
        s = ( fortran ? array::make_shape(nproma, nlev, nblk) : array::make_shape(nblk, nlev, nproma) );
    }
    else if (nvar > 0) {
        s = ( fortran ? array::make_shape(nproma, nvar, nblk) : array::make_shape(nblk, nvar, nproma) );
    }
    else {
        s = ( fortran ? array::make_shape(nproma, nblk) : array::make_shape(nblk, nproma) );
    }
    std::string name;
    params.get("name", name);
    Log::debug() << "Creating IFS " << datatype.str() << " field: " << name << "[nblk=" << nblk << "][nvar=" << nvar
                 << "][nlev=" << nlev << "][nproma=" << nproma << "]\n";

    return FieldImpl::create(name, datatype, s);
}

namespace {
static FieldCreatorBuilder<FieldCreatorIFS> __IFS("IFS");
}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
