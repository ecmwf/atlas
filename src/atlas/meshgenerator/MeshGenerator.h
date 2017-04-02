/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>
#include <string>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/config/Parametrisation.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/Config.h"

namespace eckit { class MD5; }

namespace atlas {
namespace mesh {
    class Mesh;
} }

namespace atlas {
namespace grid {
    class Distribution;
} }

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorImpl : public eckit::Owned {

public:

    MeshGeneratorImpl();

    virtual ~MeshGeneratorImpl();

    virtual void hash(eckit::MD5&) const = 0;

    virtual void generate( const grid::Grid&, const grid::Distribution&, mesh::Mesh& ) const =0;
    virtual void generate( const grid::Grid&, mesh::Mesh& ) const =0;

    mesh::Mesh generate( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh generate( const grid::Grid& ) const;

    mesh::Mesh operator()( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh operator()( const grid::Grid& ) const;

protected:

    void generate_global_element_numbering( mesh::Mesh& mesh ) const;
    void set_projection( mesh::Mesh&, const grid::Projection& ) const;
};

//----------------------------------------------------------------------------------------------------------------------
class MeshGenerator {

public:

  using meshgenerator_t = MeshGeneratorImpl;
  typedef atlas::util::Config Parameters;
  
private:
  
  eckit::SharedPtr< const meshgenerator_t > meshgenerator_;

public:

    MeshGenerator();
    MeshGenerator( const meshgenerator_t* );
    MeshGenerator( const MeshGenerator& );
    MeshGenerator(const std::string &, const eckit::Parametrisation & = util::NoConfig());

    void hash(eckit::MD5&) const;

    mesh::Mesh generate( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh generate( const grid::Grid& ) const;

    mesh::Mesh operator()( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh operator()( const grid::Grid& ) const;
    
    const meshgenerator_t* get() const { return meshgenerator_.get(); }

};

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorFactory {
public:

    /*!
     * \brief build MeshGenerator with factory key, and default options
     * \return mesh generator
     */
    static const MeshGenerator::meshgenerator_t* build(const std::string&);

    /*!
     * \brief build MeshGenerator with factory key inside parametrisation,
     * and options specified in parametrisation as well
     * \return mesh generator
     */
    static const MeshGenerator::meshgenerator_t* build(const std::string&, const eckit::Parametrisation&);

    /*!
     * \brief list all registered mesh generators
     */
    static void list(std::ostream &);

private:

    std::string name_;
    virtual const MeshGenerator::meshgenerator_t* make() = 0 ;
    virtual const MeshGenerator::meshgenerator_t* make(const eckit::Parametrisation&) = 0 ;

protected:

    MeshGeneratorFactory(const std::string&);
    virtual ~MeshGeneratorFactory();

};

//----------------------------------------------------------------------------------------------------------------------

template<class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
  virtual const MeshGenerator::meshgenerator_t* make() {
      return new T();
  }
  virtual const MeshGenerator::meshgenerator_t* make(const eckit::Parametrisation& param) {
        return new T(param);
    }
  public:
    MeshGeneratorBuilder(const std::string& name) : MeshGeneratorFactory(name) {}
};

//----------------------------------------------------------------------------------------------------------------------

extern "C" {
void atlas__MeshGenerator__delete(MeshGenerator::meshgenerator_t* This);
const MeshGenerator::meshgenerator_t* atlas__MeshGenerator__create_noconfig(const char* name);
const MeshGenerator::meshgenerator_t* atlas__MeshGenerator__create(const char* name, const eckit::Parametrisation* params);
mesh::Mesh::Implementation* atlas__MeshGenerator__generate__grid_griddist(const MeshGenerator::meshgenerator_t* This, const grid::Grid::grid_t* grid, const grid::Distribution::impl_t* distribution);
mesh::Mesh::Implementation* atlas__MeshGenerator__generate__grid(const MeshGenerator::meshgenerator_t* This, const grid::Grid::grid_t* grid);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerators
} // namespace atlas
