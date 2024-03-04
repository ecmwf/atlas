/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/nonlinear/Missing.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {



static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1(MissingIfAllMissing::static_type());
static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2(MissingIfAnyMissing::static_type());
static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3(MissingIfHeaviestMissing::static_type());

// Deprecated factory entries with "-real32" and "-real64" suffix for backwards compatibility.
static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1_real32(MissingIfAllMissing::static_type()+"-real32");
static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2_real32(MissingIfAnyMissing::static_type()+"-real32");
static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3_real32(MissingIfHeaviestMissing::static_type()+"-real32");
static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1_real64(MissingIfAllMissing::static_type()+"-real64");
static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2_real64(MissingIfAnyMissing::static_type()+"-real64");
static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3_real64(MissingIfHeaviestMissing::static_type()+"-real64");

namespace {
struct force_link {
    template <typename M>
    void load_builder() {
        NonLinearFactoryBuilder<M>("tmp");
    }
    force_link() {
        load_builder<MissingIfAllMissing>();
        load_builder<MissingIfAnyMissing>();
        load_builder<MissingIfHeaviestMissing>();
    }
};
}  // namespace
void force_link_missing() {
    static force_link static_linking;
}


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
