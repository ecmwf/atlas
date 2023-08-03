/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   ElementType.h
/// @author Willem Deconinck
/// @date   October 2015

#pragma once

#include <string>

#include "atlas/library/config.h"
#include "atlas/util/Object.h"

namespace atlas {
namespace mesh {

/**
 * \brief ElementType class (abstract) that provides access to geometric
 * information of an element
 */
class ElementType : public util::Object {
public:  // methods
    static ElementType* create(const std::string&);

    //-- Constructors

    ElementType();

    ~ElementType() = 0;

    virtual const std::string& name() const = 0;

    virtual size_t dimensionality() const = 0;

    virtual size_t nb_vertices() const = 0;

    virtual size_t nb_edges() const = 0;

    virtual size_t nb_faces() const = 0;

    virtual size_t nb_nodes() const = 0;

    virtual bool parametric() const = 0;

    virtual bool simplex() const = 0;

};

extern "C"
{
ElementType* atlas__mesh__ElementType__create(const char* name);
void atlas__mesh__ElementType__delete(ElementType* This);
idx_t atlas__mesh__ElementType__nb_nodes(const ElementType* This);
idx_t atlas__mesh__ElementType__nb_edges(const ElementType* This);
int atlas__mesh__ElementType__parametric(const ElementType* This);
const char* atlas__mesh__ElementType__name(const ElementType* This);
}

//------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

// Following includes are added for backwards compatibility with temporary elementtypes.
// They contain deprecation warnings.
// When deprecation is final, these lines can be removed.
#include "atlas/mesh/elementtypes/Quadrilateral.h"
#include "atlas/mesh/elementtypes/Triangle.h"
#include "atlas/mesh/elementtypes/Line.h"
