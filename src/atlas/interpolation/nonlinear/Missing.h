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


class Missing : public NonLinear {
public:
    bool applicable(const Field& f) const override;

    using NonLinear::execute;
    bool execute(Matrix& W, const Field& field) const override {
        return do_execute(W, field, nullptr);
    }
    bool execute(Matrix& W, const Field& field, RowIndices& missing_rows) const {
        return do_execute(W, field, &missing_rows);
    }
    bool execute(Matrix& W, const array::Array& array , const Config& config) const override {
        return do_execute(W, array, config, nullptr);
    }
    bool execute(Matrix& W, const array::Array& array , const Config& config, RowIndices& missing_rows) const {
        return do_execute(W, array, config, &missing_rows);
    }
private:
    bool do_execute(Matrix& W, const Field& field, RowIndices* missing_rows) const;
    virtual bool do_execute(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const = 0;
};


class MissingIfAllMissing : public Missing {
public:
    static std::string static_type() { return "missing-if-all-missing"; }
    using Missing::execute;
    bool do_execute(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const override;
private:
    template<typename T>
    bool executeT(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const;
};


class MissingIfAnyMissing : public Missing {
public:
    static std::string static_type() { return "missing-if-any-missing"; }
    using Missing::execute;
    bool do_execute(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const override;
private:
    template<typename T>
    bool executeT(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const;
};


class MissingIfHeaviestMissing : public Missing {
public:
    static std::string static_type() { return "missing-if-heaviest-missing"; }
    using Missing::execute;
    bool do_execute(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const override;
private:
    template<typename T>
    bool executeT(Matrix& W, const array::Array&, const Config&, RowIndices* missing_rows) const;
};

}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
