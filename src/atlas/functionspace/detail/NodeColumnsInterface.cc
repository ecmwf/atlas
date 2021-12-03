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

#include "NodeColumnsInterface.h"
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
const NodeColumns* atlas__NodesFunctionSpace__new(Mesh::Implementation* mesh, const eckit::Configuration* config) {
    ATLAS_ASSERT(mesh);
    Mesh m(mesh);
    return new NodeColumns(m, *config);
}

void atlas__NodesFunctionSpace__delete(NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    delete (This);
}

int atlas__NodesFunctionSpace__nb_nodes(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return This->nb_nodes();
}

const Mesh::Implementation* atlas__NodesFunctionSpace__mesh(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return This->mesh().get();
}

mesh::Nodes* atlas__NodesFunctionSpace__nodes(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return &This->nodes();
}

void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodeColumns* This, field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    FieldSet f(fieldset);
    This->haloExchange(f);
}

void atlas__NodesFunctionSpace__halo_exchange_field(const NodeColumns* This, field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    Field f(field);
    This->haloExchange(f);
}

const parallel::HaloExchange* atlas__NodesFunctionSpace__get_halo_exchange(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return &This->halo_exchange();
}

void atlas__NodesFunctionSpace__gather_fieldset(const NodeColumns* This, const field::FieldSetImpl* local,
                                                field::FieldSetImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet l(local);
    FieldSet g(global);
    This->gather(l, g);
}

void atlas__NodesFunctionSpace__gather_field(const NodeColumns* This, const field::FieldImpl* local,
                                             field::FieldImpl* global) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field l(local);
    Field g(global);
    This->gather(l, g);
}

const parallel::GatherScatter* atlas__NodesFunctionSpace__get_gather(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return &This->gather();
}

const parallel::GatherScatter* atlas__NodesFunctionSpace__get_scatter(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return &This->scatter();
}

void atlas__NodesFunctionSpace__scatter_fieldset(const NodeColumns* This, const field::FieldSetImpl* global,
                                                 field::FieldSetImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_FieldSet");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_FieldSet");
    const FieldSet g(global);
    FieldSet l(local);
    This->scatter(g, l);
}

void atlas__NodesFunctionSpace__scatter_field(const NodeColumns* This, const field::FieldImpl* global,
                                              field::FieldImpl* local) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(local != nullptr, "Cannot access uninitialised local atlas_Field");
    ATLAS_ASSERT(global != nullptr, "Cannot access uninitialised global atlas_Field");
    const Field g(global);
    Field l(local);
    This->scatter(g, l);
}

const parallel::Checksum* atlas__NodesFunctionSpace__get_checksum(const NodeColumns* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    return &This->checksum();
}

void atlas__NodesFunctionSpace__checksum_fieldset(const NodeColumns* This, const field::FieldSetImpl* fieldset,
                                                  char*& checksum, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    std::string checksum_str(This->checksum(fieldset));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

void atlas__NodesFunctionSpace__checksum_field(const NodeColumns* This, const field::FieldImpl* field, char*& checksum,
                                               int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::string checksum_str(This->checksum(field));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

void atlas__NodesFunctionSpace__sum_double(const NodeColumns* This, const field::FieldImpl* field, double& sum,
                                           int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_float(const NodeColumns* This, const field::FieldImpl* field, float& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_long(const NodeColumns* This, const field::FieldImpl* field, long& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_int(const NodeColumns* This, const field::FieldImpl* field, int& sum, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->sum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                               int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__sum_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                              int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__sum_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& sum,
                                             int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__sum_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& sum,
                                            int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__oisum_double(const NodeColumns* This, const field::FieldImpl* field, double& sum,
                                             int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->orderIndependentSum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_float(const NodeColumns* This, const field::FieldImpl* field, float& sum,
                                            int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->orderIndependentSum(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_long(const NodeColumns* This, const field::FieldImpl* field, long& sum, int& N) {
    atlas__NodesFunctionSpace__sum_long(This, field, sum, N);
}

void atlas__NodesFunctionSpace__oisum_int(const NodeColumns* This, const field::FieldImpl* field, int& sum, int& N) {
    atlas__NodesFunctionSpace__sum_int(This, field, sum, N);
}

void atlas__NodesFunctionSpace__oisum_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                                 int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__oisum_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                                int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__oisum_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& sum,
                                              int& size, int& N) {
    atlas__NodesFunctionSpace__sum_arr_int(This, field, sum, size, N);
}

void atlas__NodesFunctionSpace__oisum_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& sum,
                                               int& size, int& N) {
    atlas__NodesFunctionSpace__sum_arr_long(This, field, sum, size, N);
}


void atlas__NodesFunctionSpace__min_double(const NodeColumns* This, const field::FieldImpl* field, double& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__NodesFunctionSpace__min_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__NodesFunctionSpace__min_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__NodesFunctionSpace__min_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->minimum(field, minimum);
}

void atlas__NodesFunctionSpace__max_double(const NodeColumns* This, const field::FieldImpl* field, double& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__NodesFunctionSpace__max_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__NodesFunctionSpace__max_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__NodesFunctionSpace__max_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    This->maximum(field, maximum);
}

void atlas__NodesFunctionSpace__min_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& minimum,
                                               int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__NodesFunctionSpace__min_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& minimum,
                                              int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__NodesFunctionSpace__min_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                             int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    };
}

void atlas__NodesFunctionSpace__min_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                            int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> minvec;
    This->minimum(field, minvec);
    size    = minvec.size();
    minimum = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        minimum[j] = minvec[j];
    }
}

void atlas__NodesFunctionSpace__max_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& maximum,
                                               int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<double> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new double[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__NodesFunctionSpace__max_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& maximum,
                                              int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<float> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new float[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__NodesFunctionSpace__max_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                             int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<long> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new long[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__NodesFunctionSpace__max_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                            int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    std::vector<int> maxvec;
    This->maximum(field, maxvec);
    size    = maxvec.size();
    maximum = new int[size];
    for (idx_t j = 0; j < (idx_t)size; ++j) {
        maximum[j] = maxvec[j];
    }
}

void atlas__NodesFunctionSpace__minloc_double(const NodeColumns* This, const field::FieldImpl* field, double& minimum,
                                              long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                             long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                            long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                           long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->minimumAndLocation(field, minimum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_double(const NodeColumns* This, const field::FieldImpl* field, double& maximum,
                                              long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                             long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                            long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                           long& glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    This->maximumAndLocation(field, maximum, gidx);
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                  double*& minimum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloc_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                 float*& minimum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloc_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                                long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloc_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                               long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                  double*& maximum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloc_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                 float*& maximum, long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloc_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                                long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloc_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                               long*& glb_idx, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__mean_double(const NodeColumns* This, const field::FieldImpl* field, double& mean,
                                            int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->mean(field, mean, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_float(const NodeColumns* This, const field::FieldImpl* field, float& mean,
                                           int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->mean(field, mean, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_long(const NodeColumns* This, const field::FieldImpl* field, long& mean, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_int(const NodeColumns* This, const field::FieldImpl* field, int& mean, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& mean,
                                                int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__mean_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& mean,
                                               int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__mean_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& mean,
                                              int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& mean,
                                             int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                       double& mean, double& stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->meanAndStandardDeviation(field, mean, stddev, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodeColumns* This, const field::FieldImpl* field,
                                                      float& mean, float& stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    This->meanAndStandardDeviation(field, mean, stddev, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodeColumns* This, const field::FieldImpl* field, long& mean,
                                                     long& stddev, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodeColumns* This, const field::FieldImpl* field, int& mean,
                                                    int& stddev, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                           double*& mean, double*& stddev, int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                          float*& mean, float*& stddev, int& size, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                         long*& mean, long*& stddev, int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodeColumns* This, const field::FieldImpl* field,
                                                        int*& mean, int*& stddev, int& size, int& N) {
    ATLAS_NOTIMPLEMENTED;
}

void atlas__NodesFunctionSpace__minloclev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                 double& minimum, long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                                long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                               long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                              long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->minimumAndLocation(field, minimum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                 double& maximum, long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                                long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                               long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                              long& glb_idx, int& level) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    gidx_t gidx;
    idx_t lev;
    This->maximumAndLocation(field, maximum, gidx, lev);
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                     double*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloclev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                    float*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloclev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                   long*& minimum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__minloclev_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                                  long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloclev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                     double*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloclev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                    float*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloclev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                   long*& maximum, long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__maxloclev_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                                  long*& glb_idx, int*& level, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
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

void atlas__NodesFunctionSpace__sum_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* column, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(column != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    Field sum(column);
    This->sumPerLevel(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                field::FieldImpl* column, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(column != nullptr, "Cannot access uninitialised atlas_Field");
    idx_t idx_t_N;
    Field sum(column);
    This->orderIndependentSumPerLevel(field, sum, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__min_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* min) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(min != nullptr, "Cannot access uninitialised min atlas_Field");
    Field fmin(min);
    This->minimumPerLevel(field, fmin);
}

void atlas__NodesFunctionSpace__max_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* max) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(max != nullptr, "Cannot access uninitialised max atlas_Field");
    Field fmax(max);
    This->maximumPerLevel(field, fmax);
}

void atlas__NodesFunctionSpace__minloc_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* min, field::FieldImpl* glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(min != nullptr, "Cannot access uninitialised min atlas_Field");
    ATLAS_ASSERT(glb_idx != nullptr, "Cannot access uninitialised glb_idx atlas_Field");
    Field fmin(min);
    Field fglb_idx(glb_idx);
    This->minimumAndLocationPerLevel(field, fmin, fglb_idx);
}

void atlas__NodesFunctionSpace__maxloc_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* max, field::FieldImpl* glb_idx) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(max != nullptr, "Cannot access uninitialised max atlas_Field");
    ATLAS_ASSERT(glb_idx != nullptr, "Cannot access uninitialised glb_idx atlas_Field");
    Field fmax(max);
    Field fglb_idx(glb_idx);
    This->maximumAndLocationPerLevel(field, fmax, fglb_idx);
}

void atlas__NodesFunctionSpace__mean_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* mean, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(mean != nullptr, "Cannot access uninitialised mean atlas_Field");
    idx_t idx_t_N;
    Field fmean(mean);
    This->meanPerLevel(field, fmean, idx_t_N);
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                          field::FieldImpl* mean, field::FieldImpl* stddev, int& N) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_functionspace_NodeColumns");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(mean != nullptr, "Cannot access uninitialised mean atlas_Field");
    ATLAS_ASSERT(stddev);
    idx_t idx_t_N;
    Field fmean(mean);
    Field fstddev(stddev);
    This->meanAndStandardDeviationPerLevel(field, fmean, fstddev, idx_t_N);
    N = idx_t_N;
}
}

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
