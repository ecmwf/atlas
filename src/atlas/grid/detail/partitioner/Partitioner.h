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

#include <string>

#include "atlas/util/Object.h"
#include "atlas/util/Config.h"
#include "atlas/library/config.h"

namespace atlas {
class Grid;
class Mesh;
class FunctionSpace;
namespace grid {
class Distribution;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class Partitioner : public util::Object {
public:
    using Grid = atlas::Grid;

public:
    Partitioner();
    Partitioner( const idx_t nb_partitions );
    virtual ~Partitioner();

    virtual void partition( const Grid& grid, int part[] ) const = 0;

    Distribution partition( const Grid& grid ) const;

    idx_t nb_partitions() const;

    virtual std::string type() const = 0;

private:
    idx_t nb_partitions_;
};

// ------------------------------------------------------------------

class PartitionerFactory {
public:
    using Grid = Partitioner::Grid;

public:
    /*!
   * \brief build Partitioner with factory key, constructor arguments
   * \return Partitioner
   */
    static Partitioner* build( const std::string& );
    static Partitioner* build( const std::string&, const idx_t nb_partitions );
    static Partitioner* build( const std::string&, const idx_t nb_partitions, const eckit::Parametrisation & );

    /*!
   * \brief list all registered partioner builders
   */
    static void list( std::ostream& );
    static bool has( const std::string& name );

private:
    std::string name_;
    virtual Partitioner* make()                                                            = 0;
    virtual Partitioner* make( const idx_t nb_partitions )                                 = 0;
    virtual Partitioner* make( const idx_t nb_partitions, const eckit::Parametrisation & ) = 0;

protected:
    PartitionerFactory( const std::string& );
    virtual ~PartitionerFactory();
};

// ------------------------------------------------------------------

template <class T>
class PartitionerBuilder : public PartitionerFactory {
    virtual Partitioner* make() { return new T(); }

    virtual Partitioner* make( const idx_t nb_partitions ) { return new T( nb_partitions ); }
    virtual Partitioner* make( const idx_t nb_partitions, const eckit::Parametrisation & config ) { return new T( nb_partitions, config ); }

public:
    PartitionerBuilder( const std::string& name ) : PartitionerFactory( name ) {}
};

// ------------------------------------------------------------------

class MatchingPartitionerFactory {
public:
    static Partitioner* build( const std::string& type, const Mesh& partitioned );
    static Partitioner* build( const std::string& type, const FunctionSpace& partitioned );
};

// ------------------------------------------------------------------

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
