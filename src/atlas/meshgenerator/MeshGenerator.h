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
#include "eckit/config/Parametrisation.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Distribution.h"
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

class MeshGenerator : public eckit::Owned {

public:

  typedef eckit::SharedPtr<MeshGenerator> Ptr;
  typedef atlas::util::Config Parameters;

  static MeshGenerator* create(const std::string &, const eckit::Parametrisation & = Parameters());

public:

    MeshGenerator();

    virtual ~MeshGenerator();

    virtual void hash(eckit::MD5&) const = 0;

    virtual void generate( const grid::Grid&, const grid::Distribution&, mesh::Mesh& ) const =0;
    virtual void generate( const grid::Grid&, mesh::Mesh& ) const =0;

    mesh::Mesh* generate( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh* generate( const grid::Grid& ) const;

    mesh::Mesh* operator()( const grid::Grid&, const grid::Distribution& ) const;
    mesh::Mesh* operator()( const grid::Grid& ) const;

protected:

    void generate_global_element_numbering( mesh::Mesh& mesh ) const;
};

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorFactory {
public:

    /*!
     * \brief build MeshGenerator with factory key, and default options
     * \return mesh generator
     */
    static MeshGenerator* build(const std::string&);

    /*!
     * \brief build MeshGenerator with factory key inside parametrisation,
     * and options specified in parametrisation as well
     * \return mesh generator
     */
    static MeshGenerator* build(const std::string&, const eckit::Parametrisation&);

    /*!
     * \brief list all registered mesh generators
     */
    static void list(std::ostream &);

private:

    std::string name_;
    virtual MeshGenerator* make() = 0 ;
    virtual MeshGenerator* make(const eckit::Parametrisation&) = 0 ;

protected:

    MeshGeneratorFactory(const std::string&);
    virtual ~MeshGeneratorFactory();

};

//----------------------------------------------------------------------------------------------------------------------

template<class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
  virtual MeshGenerator* make() {
      return new T();
  }
  virtual MeshGenerator* make(const eckit::Parametrisation& param) {
        return new T(param);
    }
  public:
    MeshGeneratorBuilder(const std::string& name) : MeshGeneratorFactory(name) {}
};

//----------------------------------------------------------------------------------------------------------------------

extern "C" {
void atlas__MeshGenerator__delete(MeshGenerator* This);
MeshGenerator* atlas__MeshGenerator__create_noconfig(const char* name);
MeshGenerator* atlas__MeshGenerator__create(const char* name, const eckit::Parametrisation* params);
mesh::Mesh* atlas__MeshGenerator__generate__grid_griddist(const MeshGenerator* This, const grid::Grid::grid_t* grid, const grid::Distribution::impl_t* distribution);
mesh::Mesh* atlas__MeshGenerator__generate__grid(const MeshGenerator* This, const grid::Grid::grid_t* grid);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerators
} // namespace atlas
