/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library/config.h"
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif

#include "atlas/functionspace/Spectral.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/ifs/VorDivToUVIFS.h"

using atlas::FunctionSpace;
using atlas::functionspace::Spectral;

namespace atlas {
namespace trans {

namespace {
static VorDivToUVBuilder<VorDivToUVIFS> builder_ifs("ifs");  // Deprecated
static VorDivToUVBuilder<VorDivToUVIFS> builder_ectrans("ectrans");
}  // namespace

namespace {
void trans_check(const int code, const char* msg, const eckit::CodeLocation& location) {
    if (code != TRANS_SUCCESS) {
        std::stringstream errmsg;
        errmsg << "atlas::trans ERROR: " << msg << " failed: \n";
        errmsg << ::trans_error_msg(code);
        throw_Exception(errmsg.str(), location);
    }
}
#define TRANS_CHECK(CALL) trans_check(CALL, #CALL, Here())

}  // namespace

void VorDivToUVIFS::execute(const int nb_coeff, const int nb_fields, const double vorticity[],
                            const double divergence[], double U[], double V[], const eckit::Configuration&) const {
    struct ::VorDivToUV_t vordiv_to_UV = new_vordiv_to_UV();
    vordiv_to_UV.rspvor                = vorticity;
    vordiv_to_UV.rspdiv                = divergence;
    vordiv_to_UV.rspu                  = U;
    vordiv_to_UV.rspv                  = V;
    vordiv_to_UV.nfld                  = nb_fields;
    vordiv_to_UV.ncoeff                = nb_coeff;
    vordiv_to_UV.nsmax                 = truncation_;
    TRANS_CHECK(::trans_vordiv_to_UV(&vordiv_to_UV));
}

VorDivToUVIFS::VorDivToUVIFS(const int truncation, const eckit::Configuration&): truncation_(truncation) {}

VorDivToUVIFS::VorDivToUVIFS(const FunctionSpace& fs, const eckit::Configuration&):
    truncation_(Spectral(fs).truncation()) {}

VorDivToUVIFS::~VorDivToUVIFS() = default;

}  // namespace trans
}  // namespace atlas
