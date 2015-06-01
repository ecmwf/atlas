/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_meshgen_MeshGenerator_h
#define atlas_meshgen_MeshGenerator_h

#include <iosfwd>
#include <string>

#include "eckit/memory/NonCopyable.h"

namespace atlas {

class Mesh;
class Grid;
class GridDistribution;

namespace meshgen {

//------------------------------------------------------------------------------------------------------

class MeshGenerator : private eckit::NonCopyable {

public:

    MeshGenerator();

    virtual ~MeshGenerator();

    virtual void generate( const Grid&, const GridDistribution&, Mesh& ) const =0;
    virtual void generate( const Grid&, Mesh& ) const =0;

    Mesh* operator()( const Grid&, const GridDistribution& ) const;
    Mesh* operator()( const Grid& ) const;

};

class MeshGeneratorFactory {
    std::string name_;
    virtual MeshGenerator *make() = 0;

  protected:

    MeshGeneratorFactory(const std::string &);
    virtual ~MeshGeneratorFactory();

  public:

    static void list(std::ostream &);
    static MeshGenerator *build(const std::string &);

};

template< class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
    virtual MeshGenerator *make() {
        return new T();
    }
  public:
    MeshGeneratorBuilder(const std::string &name) : MeshGeneratorFactory(name) {}
};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif
