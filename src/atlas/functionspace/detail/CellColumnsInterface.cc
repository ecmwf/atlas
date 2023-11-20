/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "CellColumnsInterface.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace functionspace {
namespace detail {

using atlas::FieldSet;
using atlas::field::FieldImpl;
using atlas::field::FieldSetImpl;

// ----------------------------------------------------------------------

extern "C" {
const CellColumns* atlas__CellsFunctionSpace__new(Mesh::Implementation* mesh, const eckit::Configuration* config) {
    ATLAS_ASSERT(mesh);
    Mesh m(mesh);
    return new CellColumns(m, *config);
}

void atlas__CellsFunctionSpace__delete(CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    delete (This);
}

int atlas__CellsFunctionSpace__nb_cells(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return This->nb_cells();
}

const Mesh::Implementation* atlas__CellsFunctionSpace__mesh(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return This->mesh().get();
}

const mesh::HybridElements* atlas__CellsFunctionSpace__cells(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->cells();
}

void atlas__CellsFunctionSpace__halo_exchange_fieldset(const CellColumns* This, field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    FieldSet f(fieldset);
    This->haloExchange(f);
}

void atlas__CellsFunctionSpace__halo_exchange_field(const CellColumns* This, field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    Field f(field);
    This->haloExchange(f);
}

const parallel::HaloExchange* atlas__CellsFunctionSpace__get_halo_exchange(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->halo_exchange();
}

void atlas__CellsFunctionSpace__gather_fieldset(const CellColumns* This, const field::FieldSetImpl* local,
                                                field::FieldSetImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet l(local);
    FieldSet g(global);
    This->gather(l, g);
}

void atlas__CellsFunctionSpace__gather_field(const CellColumns* This, const field::FieldImpl* local,
                                             field::FieldImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field l(local);
    Field g(global);
    This->gather(l, g);
}

const parallel::GatherScatter* atlas__CellsFunctionSpace__get_gather(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->gather();
}

const parallel::GatherScatter* atlas__CellsFunctionSpace__get_scatter(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->scatter();
}

void atlas__CellsFunctionSpace__scatter_fieldset(const CellColumns* This, const field::FieldSetImpl* global,
                                                 field::FieldSetImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet g(global);
    FieldSet l(local);
    This->scatter(g, l);
}

void atlas__CellsFunctionSpace__scatter_field(const CellColumns* This, const field::FieldImpl* global,
                                              field::FieldImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field g(global);
    Field l(local);
    This->scatter(g, l);
}

const parallel::Checksum* atlas__CellsFunctionSpace__get_checksum(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    return &This->checksum();
}

void atlas__CellsFunctionSpace__checksum_fieldset(const CellColumns* This, const field::FieldSetImpl* fieldset,
                                                  char*& checksum, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    std::string checksum_str(This->checksum(fieldset));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

void atlas__CellsFunctionSpace__checksum_field(const CellColumns* This, const field::FieldImpl* field, char*& checksum,
                                               int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::string checksum_str(This->checksum(field));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

/*
void atlas__CellsFunctionSpace__sum_double(const CellColumns* This, const field::FieldImpl* field, double& sum,
                                           int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_float(const CellColumns* This, const field::FieldImpl* field, float& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_long(const CellColumns* This, const field::FieldImpl* field, long& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_int(const CellColumns* This, const field::FieldImpl* field, int& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_arr_double(const CellColumns* This, const field::FieldImpl* field, double*& sum,
                                               int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<double> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_arr_float(const CellColumns* This, const field::FieldImpl* field, float*& sum,
                                              int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<float> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& sum,
                                             int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<long> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__sum_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& sum,
                                            int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<int> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_double(const CellColumns* This, const field::FieldImpl* field, double& sum,
                                             int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->orderIndependentSum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_float(const CellColumns* This, const field::FieldImpl* field, float& sum,
                                            int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->orderIndependentSum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_long(const CellColumns* This, const field::FieldImpl* field, long& sum, int& N) {
    atlas__CellsFunctionSpace__sum_long(This, field, sum, N);
}

void atlas__CellsFunctionSpace__oisum_int(const CellColumns* This, const field::FieldImpl* field, int& sum, int& N) {
    atlas__CellsFunctionSpace__sum_int(This, field, sum, N);
}

void atlas__CellsFunctionSpace__oisum_arr_double(const CellColumns* This, const field::FieldImpl* field, double*& sum,
                                                 int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<double> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_arr_float(const CellColumns* This, const field::FieldImpl* field, float*& sum,
                                                int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<float> sumvec;
    This->orderIndependentSum(field, sumvec, idx_t_N);
    size = sumvec.size();
    sum  = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        sum[j] = sumvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& sum,
                                              int& size, int& N) {
    atlas__CellsFunctionSpace__sum_arr_int(This, field, sum, size, N);
}

void atlas__CellsFunctionSpace__oisum_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& sum,
                                               int& size, int& N) {
    atlas__CellsFunctionSpace__sum_arr_long(This, field, sum, size, N);
}


void atlas__CellsFunctionSpace__min_double(const CellColumns* This, const field::FieldImpl* field, double& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__CellsFunctionSpace__min_float(const CellColumns* This, const field::FieldImpl* field, float& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__CellsFunctionSpace__min_long(const CellColumns* This, const field::FieldImpl* field, long& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__CellsFunctionSpace__min_int(const CellColumns* This, const field::FieldImpl* field, int& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__CellsFunctionSpace__max_double(const CellColumns* This, const field::FieldImpl* field, double& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__CellsFunctionSpace__max_float(const CellColumns* This, const field::FieldImpl* field, float& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__CellsFunctionSpace__max_long(const CellColumns* This, const field::FieldImpl* field, long& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__CellsFunctionSpace__max_int(const CellColumns* This, const field::FieldImpl* field, int& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__CellsFunctionSpace__min_arr_double(const CellColumns* This, const field::FieldImpl* field, double*& minimum,
                                               int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__CellsFunctionSpace__min_arr_float(const CellColumns* This, const field::FieldImpl* field, float*& minimum,
                                              int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__CellsFunctionSpace__min_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& minimum,
                                             int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    };
}

void atlas__CellsFunctionSpace__min_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& minimum,
                                            int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__CellsFunctionSpace__max_arr_double(const CellColumns* This, const field::FieldImpl* field, double*& maximum,
                                               int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__CellsFunctionSpace__max_arr_float(const CellColumns* This, const field::FieldImpl* field, float*& maximum,
                                              int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__CellsFunctionSpace__max_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& maximum,
                                             int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__CellsFunctionSpace__max_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& maximum,
                                            int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__CellsFunctionSpace__minloc_double(const CellColumns* This, const field::FieldImpl* field, double& minimum,
                                              long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__minloc_float(const CellColumns* This, const field::FieldImpl* field, float& minimum,
                                             long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__minloc_long(const CellColumns* This, const field::FieldImpl* field, long& minimum,
                                            long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__minloc_int(const CellColumns* This, const field::FieldImpl* field, int& minimum,
                                           long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__maxloc_double(const CellColumns* This, const field::FieldImpl* field, double& maximum,
                                              long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__maxloc_float(const CellColumns* This, const field::FieldImpl* field, float& maximum,
                                             long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__maxloc_long(const CellColumns* This, const field::FieldImpl* field, long& maximum,
                                            long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__maxloc_int(const CellColumns* This, const field::FieldImpl* field, int& maximum,
                                           long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__CellsFunctionSpace__minloc_arr_double(const CellColumns* This, const field::FieldImpl* field,
                                                  double*& minimum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> minvec;
    std::vector<gidx_t> gidxvec;
    This->minimumAndLocation(field, minvec, gidxvec);
    size    = minvec.size();
    minimum = new double[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__minloc_arr_float(const CellColumns* This, const field::FieldImpl* field,
                                                 float*& minimum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> minvec;
    std::vector<gidx_t> gidxvec;
    This->minimumAndLocation(field, minvec, gidxvec);
    size    = minvec.size();
    minimum = new float[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__minloc_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& minimum,
                                                long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> minvec;
    std::vector<gidx_t> gidxvec;
    This->minimumAndLocation(field, minvec, gidxvec);
    size    = minvec.size();
    minimum = new long[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__minloc_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& minimum,
                                               long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> minvec;
    std::vector<gidx_t> gidxvec;
    This->minimumAndLocation(field, minvec, gidxvec);
    size    = minvec.size();
    minimum = new int[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloc_arr_double(const CellColumns* This, const field::FieldImpl* field,
                                                  double*& maximum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> maxvec;
    std::vector<gidx_t> gidxvec;
    This->maximumAndLocation(field, maxvec, gidxvec);
    size    = maxvec.size();
    maximum = new double[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloc_arr_float(const CellColumns* This, const field::FieldImpl* field,
                                                 float*& maximum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> maxvec;
    std::vector<gidx_t> gidxvec;
    This->maximumAndLocation(field, maxvec, gidxvec);
    size    = maxvec.size();
    maximum = new float[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloc_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& maximum,
                                                long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> maxvec;
    std::vector<gidx_t> gidxvec;
    This->maximumAndLocation(field, maxvec, gidxvec);
    size    = maxvec.size();
    maximum = new long[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloc_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& maximum,
                                               long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> maxvec;
    std::vector<gidx_t> gidxvec;
    This->maximumAndLocation(field, maxvec, gidxvec);
    size    = maxvec.size();
    maximum = new int[size];
    glb_idx = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
    }
}

void atlas__CellsFunctionSpace__mean_double(const CellColumns* This, const field::FieldImpl* field, double& mean,
                                            int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->mean(field, mean, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_float(const CellColumns* This, const field::FieldImpl* field, float& mean,
                                           int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->mean(field, mean, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_long(const CellColumns* This, const field::FieldImpl* field, long& mean, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_int(const CellColumns* This, const field::FieldImpl* field, int& mean, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_arr_double(const CellColumns* This, const field::FieldImpl* field, double*& mean,
                                                int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<double> meanvec;
    This->mean(field, meanvec, idx_t_N);
    size = meanvec.size();
    mean = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        mean[j] = meanvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_arr_float(const CellColumns* This, const field::FieldImpl* field, float*& mean,
                                               int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<float> meanvec;
    This->mean(field, meanvec, idx_t_N);
    size = meanvec.size();
    mean = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        mean[j] = meanvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_arr_long(const CellColumns* This, const field::FieldImpl* field, long*& mean,
                                              int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& mean,
                                             int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_and_stddev_double(const CellColumns* This, const field::FieldImpl* field,
                                                       double& mean, double& stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->meanAndStandardDeviation(field, mean, stddev, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_and_stddev_float(const CellColumns* This, const field::FieldImpl* field,
                                                      float& mean, float& stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->meanAndStandardDeviation(field, mean, stddev, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_and_stddev_long(const CellColumns* This, const field::FieldImpl* field, long& mean,
                                                     long& stddev, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_and_stddev_int(const CellColumns* This, const field::FieldImpl* field, int& mean,
                                                    int& stddev, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_and_stddev_arr_double(const CellColumns* This, const field::FieldImpl* field,
                                                           double*& mean, double*& stddev, int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<double> meanvec;
    std::vector<double> stddevvec;
    This->meanAndStandardDeviation(field, meanvec, stddevvec, idx_t_N);
    size   = meanvec.size();
    mean   = new double[size];
    stddev = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        mean[j]   = meanvec[j];
        stddev[j] = stddevvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_and_stddev_arr_float(const CellColumns* This, const field::FieldImpl* field,
                                                          float*& mean, float*& stddev, int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    std::vector<float> meanvec;
    std::vector<float> stddevvec;
    This->meanAndStandardDeviation(field, meanvec, stddevvec, idx_t_N);
    size   = meanvec.size();
    mean   = new float[size];
    stddev = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        mean[j]   = meanvec[j];
        stddev[j] = stddevvec[j];
    }
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_and_stddev_arr_long(const CellColumns* This, const field::FieldImpl* field,
                                                         long*& mean, long*& stddev, int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__mean_and_stddev_arr_int(const CellColumns* This, const field::FieldImpl* field,
                                                        int*& mean, int*& stddev, int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__CellsFunctionSpace__minloclev_double(const CellColumns* This, const field::FieldImpl* field,
                                                 double& minimum, long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__minloclev_float(const CellColumns* This, const field::FieldImpl* field, float& minimum,
                                                long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__minloclev_long(const CellColumns* This, const field::FieldImpl* field, long& minimum,
                                               long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__minloclev_int(const CellColumns* This, const field::FieldImpl* field, int& minimum,
                                              long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__maxloclev_double(const CellColumns* This, const field::FieldImpl* field,
                                                 double& maximum, long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__maxloclev_float(const CellColumns* This, const field::FieldImpl* field, float& maximum,
                                                long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__maxloclev_long(const CellColumns* This, const field::FieldImpl* field, long& maximum,
                                               long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__maxloclev_int(const CellColumns* This, const field::FieldImpl* field, int& maximum,
                                              long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__CellsFunctionSpace__minloclev_arr_double(const CellColumns* This, const field::FieldImpl* field,
                                                     double*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> minvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->minimumAndLocation(field, minvec, gidxvec, levvec);
    size    = minvec.size();
    minimum = new double[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__minloclev_arr_float(const CellColumns* This, const field::FieldImpl* field,
                                                    float*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> minvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->minimumAndLocation(field, minvec, gidxvec, levvec);
    size    = minvec.size();
    minimum = new float[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    };
}

void atlas__CellsFunctionSpace__minloclev_arr_long(const CellColumns* This, const field::FieldImpl* field,
                                                   long*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> minvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->minimumAndLocation(field, minvec, gidxvec, levvec);
    size    = minvec.size();
    minimum = new long[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__minloclev_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& minimum,
                                                  long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> minvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->minimumAndLocation(field, minvec, gidxvec, levvec);
    size    = minvec.size();
    minimum = new int[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    };
}

void atlas__CellsFunctionSpace__maxloclev_arr_double(const CellColumns* This, const field::FieldImpl* field,
                                                     double*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> maxvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->maximumAndLocation(field, maxvec, gidxvec, levvec);
    size    = maxvec.size();
    maximum = new double[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloclev_arr_float(const CellColumns* This, const field::FieldImpl* field,
                                                    float*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> maxvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->maximumAndLocation(field, maxvec, gidxvec, levvec);
    size    = maxvec.size();
    maximum = new float[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloclev_arr_long(const CellColumns* This, const field::FieldImpl* field,
                                                   long*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> maxvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->maximumAndLocation(field, maxvec, gidxvec, levvec);
    size    = maxvec.size();
    maximum = new long[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__maxloclev_arr_int(const CellColumns* This, const field::FieldImpl* field, int*& maximum,
                                                  long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> maxvec;
    std::vector<gidx_t> gidxvec;
    std::vector<idx_t> levvec;
    This->maximumAndLocation(field, maxvec, gidxvec, levvec);
    size    = maxvec.size();
    maximum = new int[size];
    glb_idx = new long[size];
    level   = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
        glb_idx[j] = gidxvec[j];
        level[j]   = levvec[j];
    }
}

void atlas__CellsFunctionSpace__sum_per_level(const CellColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* column, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(column != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    Field sum(column);
    This->sumPerLevel(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__oisum_per_level(const CellColumns* This, const field::FieldImpl* field,
                                                field::FieldImpl* column, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(column != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    Field sum(column);
    This->orderIndependentSumPerLevel(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__min_per_level(const CellColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* min) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(min != nullptr, "Cannot access uninitialised min atlas_Field");
    Field fmin(min);
    This->minimumPerLevel(field, fmin);
}

void atlas__CellsFunctionSpace__max_per_level(const CellColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* max) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(max != nullptr, "Cannot access uninitialised max atlas_Field");
    Field fmax(max);
    This->maximumPerLevel(field, fmax);
}

void atlas__CellsFunctionSpace__minloc_per_level(const CellColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* min, field::FieldImpl* glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(min != nullptr, "Cannot access uninitialised min atlas_Field");
    ATLAS_ASSERT(glb_idx != nullptr, "Cannot access uninitialised glb_idx atlas_Field");
    Field fmin(min);
    Field fglb_idx(glb_idx);
    This->minimumAndLocationPerLevel(field, fmin, fglb_idx);
}

void atlas__CellsFunctionSpace__maxloc_per_level(const CellColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* max, field::FieldImpl* glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(max != nullptr, "Cannot access uninitialised max atlas_Field");
    ATLAS_ASSERT(glb_idx != nullptr, "Cannot access uninitialised glb_idx atlas_Field");
    Field fmax(max);
    Field fglb_idx(glb_idx);
    This->maximumAndLocationPerLevel(field, fmax, fglb_idx);
}

void atlas__CellsFunctionSpace__mean_per_level(const CellColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* mean, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(mean != nullptr, "Cannot access uninitialised mean atlas_Field");
    idx_t idx_t_N;
    Field fmean(mean);
    This->meanPerLevel(field, fmean, idx_t_N);
    N = idx_t_N;
}

void atlas__CellsFunctionSpace__mean_and_stddev_per_level(const CellColumns* This, const field::FieldImpl* field,
                                                          field::FieldImpl* mean, field::FieldImpl* stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_CellColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(mean != nullptr, "Cannot access uninitialised mean atlas_Field");
    ATLAS_ASSERT(stddev);
    idx_t idx_t_N;
    Field fmean(mean);
    Field fstddev(stddev);
    This->meanAndStandardDeviationPerLevel(field, fmean, fstddev, idx_t_N);
    N = idx_t_N;
}
*/
}

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
