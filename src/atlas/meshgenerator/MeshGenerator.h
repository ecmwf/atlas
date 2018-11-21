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

#include <iosfwd>
#include <string>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/util/Config.h"

namespace eckit {
class Hash;
class Parametrisation;
}  // namespace eckit

namespace atlas {
class Mesh;
class Grid;
class Projection;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
class DistributionImpl;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
class Grid;
}  // namespace grid
}  // namespace detail
}  // namespace grid
using GridImpl = grid::detail::grid::Grid;
}  // namespace atlas
namespace atlas {
namespace mesh {
namespace detail {
class MeshImpl;
}  // namespace detail
}  // namespace mesh
}  // namespace atlas


namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorImpl : public eckit::Owned {
public:
    MeshGeneratorImpl();

    virtual ~MeshGeneratorImpl();

    virtual void hash( eckit::Hash& ) const = 0;

    virtual void generate( const Grid&, const grid::Distribution&, Mesh& ) const = 0;
    virtual void generate( const Grid&, Mesh& ) const                            = 0;

    Mesh generate( const Grid&, const grid::Distribution& ) const;
    Mesh generate( const Grid& ) const;

    Mesh operator()( const Grid&, const grid::Distribution& ) const;
    Mesh operator()( const Grid& ) const;

protected:
    void generateGlobalElementNumbering( Mesh& mesh ) const;
    void setProjection( Mesh&, const Projection& ) const;
    void setGrid( Mesh&, const Grid&, const grid::Distribution& ) const;
};

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorFactory {
public:
    /*!
   * \brief build MeshGenerator with factory key, and default options
   * \return mesh generator
   */
    static const MeshGeneratorImpl* build( const std::string& );

    /*!
   * \brief build MeshGenerator with factory key inside parametrisation,
   * and options specified in parametrisation as well
   * \return mesh generator
   */
    static const MeshGeneratorImpl* build( const std::string&, const eckit::Parametrisation& );

    /*!
   * \brief list all registered mesh generators
   */
    static void list( std::ostream& );

private:
    std::string name_;
    virtual const MeshGeneratorImpl* make()                                = 0;
    virtual const MeshGeneratorImpl* make( const eckit::Parametrisation& ) = 0;

protected:
    MeshGeneratorFactory( const std::string& );
    virtual ~MeshGeneratorFactory();
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
    virtual const MeshGeneratorImpl* make() { return new T(); }
    virtual const MeshGeneratorImpl* make( const eckit::Parametrisation& param ) { return new T( param ); }

public:
    MeshGeneratorBuilder( const std::string& name ) : MeshGeneratorFactory( name ) {}
};

//----------------------------------------------------------------------------------------------------------------------

extern "C" {
void atlas__MeshGenerator__delete( MeshGeneratorImpl* This );
const MeshGeneratorImpl* atlas__MeshGenerator__create_noconfig( const char* name );
const MeshGeneratorImpl* atlas__MeshGenerator__create( const char* name, const eckit::Parametrisation* params );
mesh::detail::MeshImpl* atlas__MeshGenerator__generate__grid_griddist( const MeshGeneratorImpl* This,
                                                                       const GridImpl* grid,
                                                                       const grid::DistributionImpl* distribution );
mesh::detail::MeshImpl* atlas__MeshGenerator__generate__grid( const MeshGeneratorImpl* This, const GridImpl* grid );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

class MeshGenerator {
public:
    using Implementation = meshgenerator::MeshGeneratorImpl;
    typedef atlas::util::Config Parameters;

private:
    eckit::SharedPtr<const Implementation> meshgenerator_;

public:
    MeshGenerator();
    MeshGenerator( const Implementation* );
    MeshGenerator( const MeshGenerator& );
    MeshGenerator( const std::string&, const eckit::Parametrisation& = util::NoConfig() );

    void hash( eckit::Hash& ) const;

    Mesh generate( const Grid&, const grid::Distribution& ) const;
    Mesh generate( const Grid& ) const;

    Mesh operator()( const Grid&, const grid::Distribution& ) const;
    Mesh operator()( const Grid& ) const;

    const Implementation* get() const { return meshgenerator_.get(); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
