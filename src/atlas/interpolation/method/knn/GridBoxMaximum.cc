/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/GridBoxMaximum.h"

#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {
MethodBuilder<GridBoxMaximum> __builder( "grid-box-maximum" );
}


void GridBoxMaximum::do_execute( const FieldSet& source, FieldSet& target ) const {
    ATLAS_ASSERT( source.size() == target.size() );

    // Matrix-based interpolation is handled by base (Method) class
    // TODO: exploit sparse/dense matrix multiplication
    for ( idx_t i = 0; i < source.size(); ++i ) {
        if ( matrixFree_ ) {
            GridBoxMaximum::do_execute( source[i], target[i] );
        }
        else {
            Method::do_execute( source[i], target[i] );
        }
    }
}


void GridBoxMaximum::do_execute( const Field& /*source*/, Field& /*target*/ ) const {
    ATLAS_NOTIMPLEMENTED;
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
