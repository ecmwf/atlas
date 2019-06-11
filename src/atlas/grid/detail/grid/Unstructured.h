/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date January 2015

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/util/Point.h"

namespace atlas {
class Mesh;
}

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Unstructured : public Grid {
public:
    class IteratorXY : public Grid::IteratorXY {
    public:
        IteratorXY( const Unstructured& grid, bool begin = true ) :
            grid_( grid ),
            size_( static_cast<idx_t>( grid_.points_->size() ) ),
            n_( begin ? 0 : size_ ) {}

        virtual bool next( PointXY& xy ) {
            if ( n_ != size_ ) {
                xy = grid_.xy( n_++ );
                return true;
            }
            else {
                return false;
            }
        }

        virtual const PointXY operator*() const { return grid_.xy( n_ ); }

        virtual const Grid::IteratorXY& operator++() {
            ++n_;
            return *this;
        }

        virtual bool operator==( const Grid::IteratorXY& other ) const {
            return n_ == static_cast<const IteratorXY&>( other ).n_;
        }

        virtual bool operator!=( const Grid::IteratorXY& other ) const {
            return n_ != static_cast<const IteratorXY&>( other ).n_;
        }

    private:
        const Unstructured& grid_;
        idx_t size_;
        idx_t n_;
    };

    class IteratorXYPredicated : public Grid::IteratorXY {
    public:
        IteratorXYPredicated( const Unstructured& grid, Grid::IteratorXY::Predicate p, bool begin = true ) :
            grid_( grid ),
            p_( p ),
            size_( static_cast<idx_t>( grid_.points_->size() ) ),
            n_( begin ? 0 : size_ ) {
            if ( begin ) {}
        }

        virtual bool next( PointXY& xy );

        virtual const PointXY operator*() const { return grid_.xy( n_ ); }

        virtual const Grid::IteratorXY& operator++() {
            do {
                ++n_;
                if ( n_ == size_ ) return *this;
            } while ( not p_( n_ ) );
            return *this;
        }

        virtual bool operator==( const Grid::IteratorXY& other ) const {
            return n_ == static_cast<const IteratorXYPredicated&>( other ).n_;
        }

        virtual bool operator!=( const Grid::IteratorXY& other ) const {
            return n_ != static_cast<const IteratorXYPredicated&>( other ).n_;
        }

    private:
        const Unstructured& grid_;
        Grid::IteratorXY::Predicate p_;
        idx_t size_;
        idx_t n_;
    };

    class IteratorLonLat : public Grid::IteratorLonLat {
    public:
        IteratorLonLat( const Unstructured& grid, bool begin = true ) :
            grid_( grid ),
            size_( static_cast<idx_t>( grid_.points_->size() ) ),
            n_( begin ? 0 : size_ ) {}

        virtual bool next( PointLonLat& lonlat ) {
            if ( n_ != size_ ) {
                lonlat = grid_.lonlat( n_++ );
                return true;
            }
            else {
                return false;
            }
        }

        virtual const PointLonLat operator*() const { return grid_.lonlat( n_ ); }

        virtual const Grid::IteratorLonLat& operator++() {
            ++n_;
            return *this;
        }

        virtual bool operator==( const Grid::IteratorLonLat& other ) const {
            return n_ == static_cast<const IteratorLonLat&>( other ).n_;
        }

        virtual bool operator!=( const Grid::IteratorLonLat& other ) const {
            return n_ != static_cast<const IteratorLonLat&>( other ).n_;
        }

    private:
        const Unstructured& grid_;
        idx_t size_;
        idx_t n_;
    };


public:  // methods
    static std::string static_type() { return "unstructured"; }
    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

    /// Constructor converting any Grid with domain to an unstructured grid
    Unstructured( const Grid&, Domain );

    /// Constructor taking a list of parameters
    Unstructured( const Config& );

    /// Constructor taking a list of points (takes ownership)
    Unstructured( std::vector<PointXY>* pts );

    /// Constructor taking a list of points (takes ownership)
    Unstructured( std::vector<PointXY>&& pts );

    /// Constructor taking a mesh
    Unstructured( const Mesh& m );

    /// Constructor from initializer list
    Unstructured( std::initializer_list<PointXY> );

    virtual ~Unstructured();

    virtual idx_t size() const;

    virtual Spec spec() const;

    PointXY xy( idx_t n ) const { return ( *points_ )[n]; }

    PointLonLat lonlat( idx_t n ) const { return projection_.lonlat( ( *points_ )[n] ); }

    virtual IteratorXY* xy_begin() const { return new IteratorXY( *this ); }
    virtual IteratorXY* xy_end() const { return new IteratorXY( *this, false ); }
    virtual IteratorLonLat* lonlat_begin() const { return new IteratorLonLat( *this ); }
    virtual IteratorLonLat* lonlat_end() const { return new IteratorLonLat( *this, false ); }
    virtual IteratorXYPredicated* xy_begin( IteratorXY::Predicate p ) const {
        return new IteratorXYPredicated( *this, p );
    }
    virtual IteratorXYPredicated* xy_end( IteratorXY::Predicate p ) const {
        return new IteratorXYPredicated( *this, p, false );
    }

private:  // methods
    virtual void print( std::ostream& ) const;

    /// Hash of the lonlat array + BoundingBox
    virtual void hash( eckit::Hash& ) const;

protected:
    /// Storage of coordinate points
    std::unique_ptr<std::vector<PointXY>> points_;

    /// Cache for the shortName
    mutable std::string shortName_;

    /// Cache for the spec since may be quite heavy to compute
    mutable std::unique_ptr<Grid::Spec> cached_spec_;
};

extern "C" {
const Unstructured* atlas__grid__Unstructured__points( const double lonlat[], int shapef[], int stridesf[] );
const Unstructured* atlas__grid__Unstructured__config( util::Config* conf );
}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
