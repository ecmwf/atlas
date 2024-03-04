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


#pragma once

#include "atlas/interpolation/nonlinear/NonLinear.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


struct Missing : NonLinear {
private:
    bool applicable(const Field& f) const override;
};


struct MissingIfAllMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const override;

    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const;

    static std::string static_type() { return "missing-if-all-missing"; }
};


struct MissingIfAnyMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const override;

    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const;

    static std::string static_type() { return "missing-if-any-missing"; }
};


struct MissingIfHeaviestMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const override;

    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const;

    static std::string static_type() { return "missing-if-heaviest-missing"; }
};

}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
