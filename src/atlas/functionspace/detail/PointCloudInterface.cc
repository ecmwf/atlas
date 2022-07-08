/*
 * (C) Copyright 200 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "PointCloudInterface.h"

#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------

extern "C" {

const detail::PointCloud* atlas__functionspace__PointCloud__new__lonlat(const Field::Implementation* lonlat) {
    return new detail::PointCloud(Field(lonlat));
}

const detail::PointCloud* atlas__functionspace__PointCloud__new__lonlat_ghost(const field::FieldImpl* lonlat,
                                                                              const field::FieldImpl* ghost) {
    return new detail::PointCloud(Field(lonlat), Field(ghost));
}


const detail::PointCloud* atlas__functionspace__PointCloud__new__grid(const Grid::Implementation* grid) {
    return new detail::PointCloud(Grid(grid));
}

const field::FieldImpl* atlas__fs__PointCloud__lonlat(const detail::PointCloud* This) {
    return This->lonlat().get();
}

idx_t atlas__fs__PointCloud__size(const detail::PointCloud* This) {
    return This->size();
}
}

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
