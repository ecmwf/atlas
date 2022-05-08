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

    //-- Accessors

    virtual const std::string& name() const = 0;
    // virtual idx_t dimensionality() const = 0;

    // virtual idx_t nb_vertices() const = 0;
    virtual idx_t nb_edges() const = 0;
    // virtual idx_t nb_faces() const = 0;

    virtual idx_t nb_nodes() const = 0;

    virtual bool parametric() const = 0;
};

namespace temporary {

class Volume : public ElementType {
public:
    enum
    {
        DIMENSIONALITY = 3
    };
};

class Face : public ElementType {
public:
    enum
    {
        DIMENSIONALITY = 2
    };
    enum
    {
        FACES = 1
    };
    virtual idx_t nb_faces() const { return FACES; }
};

class Edge : public ElementType {
public:
    enum
    {
        DIMENSIONALITY = 1
    };
    enum
    {
        FACES = 0
    };
    enum
    {
        EDGES = 1
    };
    virtual idx_t nb_faces() const { return FACES; }
    virtual idx_t nb_edges() const { return EDGES; }
};

class Vertex : public ElementType {
public:
    enum
    {
        DIMENSIONALITY = 0
    };
    enum
    {
        FACES = 0
    };
    enum
    {
        EDGES = 0
    };
    enum
    {
        VERTICES = 1
    };
    virtual idx_t nb_faces() const { return FACES; }
    virtual idx_t nb_edges() const { return EDGES; }
    virtual idx_t nb_vertices() const { return VERTICES; }
};

class Quadrilateral : public Face {
public:
    enum
    {
        EDGES = 4
    };
    enum
    {
        VERTICES = 4
    };
    enum
    {
        FACETS = EDGES
    };
    enum
    {
        RIDGES = VERTICES
    };
    virtual ~Quadrilateral() {}
    virtual bool parametric() const { return true; }
    virtual idx_t nb_vertices() const { return VERTICES; }
    virtual idx_t nb_edges() const { return EDGES; }
    virtual idx_t nb_nodes() const { return VERTICES; }
    virtual idx_t nb_facets() const { return FACETS; }
    virtual idx_t nb_ridges() const { return RIDGES; }
    virtual const std::string& name() const {
        static std::string s("Quadrilateral");
        return s;
    }
};

class Triangle : public Face {
public:
    enum
    {
        EDGES = 3
    };
    enum
    {
        VERTICES = 3
    };
    enum
    {
        FACETS = EDGES
    };
    enum
    {
        RIDGES = VERTICES
    };
    virtual ~Triangle() {}
    virtual bool parametric() const { return true; }
    virtual idx_t nb_vertices() const { return VERTICES; }
    virtual idx_t nb_edges() const { return EDGES; }
    virtual idx_t nb_nodes() const { return VERTICES; }
    virtual idx_t nb_facets() const { return FACETS; }
    virtual idx_t nb_ridges() const { return RIDGES; }
    virtual const std::string& name() const {
        static std::string s("Triangle");
        return s;
    }
};

class Line : public Edge {
public:
    enum
    {
        VERTICES = 2
    };
    enum
    {
        FACETS = VERTICES
    };
    enum
    {
        RIDGES = 0
    };
    virtual ~Line() {}
    virtual bool parametric() const { return true; }
    virtual idx_t nb_vertices() const { return VERTICES; }
    virtual idx_t nb_edges() const { return EDGES; }
    virtual idx_t nb_nodes() const { return VERTICES; }
    virtual idx_t nb_facets() const { return FACETS; }
    virtual const std::string& name() const {
        static std::string s("Line");
        return s;
    }
};

class Pentagon : public Face {
public:
    enum
    {
        EDGES = 5
    };
    enum
    {
        VERTICES = 5
    };
    enum
    {
        FACETS = EDGES
    };
    enum
    {
        RIDGES = VERTICES
    };
    virtual ~Pentagon() {}
    virtual bool parametric() const { return false; }
    virtual idx_t nb_vertices() const { return VERTICES; }
    virtual idx_t nb_edges() const { return EDGES; }
    virtual idx_t nb_nodes() const { return VERTICES; }
    virtual idx_t nb_facets() const { return FACETS; }
    virtual idx_t nb_ridges() const { return RIDGES; }
    virtual const std::string& name() const {
        static std::string s("Pentagon");
        return s;
    }
};

}  // namespace temporary

extern "C" {
ElementType* atlas__mesh__Triangle__create();
ElementType* atlas__mesh__Quadrilateral__create();
ElementType* atlas__mesh__Line__create();

void atlas__mesh__ElementType__delete(ElementType* This);
idx_t atlas__mesh__ElementType__nb_nodes(const ElementType* This);
idx_t atlas__mesh__ElementType__nb_edges(const ElementType* This);
int atlas__mesh__ElementType__parametric(const ElementType* This);
const char* atlas__mesh__ElementType__name(const ElementType* This);
}

//------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
