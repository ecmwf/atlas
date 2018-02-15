/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include "atlas/numerics/Nabla.h"

namespace atlas {
namespace numerics {
namespace fvm {
class Method;
}
}  // namespace numerics
}  // namespace atlas

namespace atlas {
class Field;
}

namespace atlas {
namespace numerics {
namespace fvm {

class Nabla : public atlas::numerics::Nabla::nabla_t {
public:
    Nabla( const atlas::numerics::Method&, const eckit::Parametrisation& );
    virtual ~Nabla();

    void gradient( const Field& scalar, Field& grad ) const;
    void divergence( const Field& vector, Field& div ) const;
    void curl( const Field& vector, Field& curl ) const;
    void laplacian( const Field& scalar, Field& laplacian ) const;

private:
    void setup();

    void gradient_of_scalar( const Field& scalar, Field& grad ) const;
    void gradient_of_vector( const Field& vector, Field& grad ) const;

private:
    fvm::Method const* fvm_;
    std::vector<size_t> pole_edges_;
};

// ------------------------------------------------------------------

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
