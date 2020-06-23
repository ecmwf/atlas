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

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

namespace atlas {
class Grid;

namespace functionspace {

//------------------------------------------------------------------------------------------------------

namespace detail {

class PointCloud : public functionspace::FunctionSpaceImpl {
public:
    PointCloud( const std::vector<PointXY>& );
    PointCloud( PointXY, const std::vector<PointXY>& );
    PointCloud( PointXYZ, const std::vector<PointXYZ>& );
    PointCloud( const Field& lonlat );
    PointCloud( const Field& lonlat, const Field& ghost );
    PointCloud( const Grid& );
    virtual ~PointCloud() override {}
    virtual std::string type() const override { return "PointCloud"; }
    virtual operator bool() const override { return true; }
    virtual size_t footprint() const override { return sizeof( *this ); }
    virtual std::string distribution() const override;
    Field lonlat() const override { return lonlat_; }
    const Field& vertical() const { return vertical_; }
    Field ghost() const override;
    virtual idx_t size() const override { return lonlat_.shape( 0 ); }

    using FunctionSpaceImpl::createField;
    virtual Field createField( const eckit::Configuration& ) const override;
    virtual Field createField( const Field&, const eckit::Configuration& ) const override;


    class IteratorXYZ {
    public:
        IteratorXYZ( const PointCloud& fs, bool begin = true );

        bool next( PointXYZ& xyz );

        const PointXYZ operator*() const;

        const IteratorXYZ& operator++() {
            ++n_;
            return *this;
        }

        bool operator==( const IteratorXYZ& other ) const { return n_ == static_cast<const IteratorXYZ&>( other ).n_; }

        virtual bool operator!=( const IteratorXYZ& other ) const {
            return n_ != static_cast<const IteratorXYZ&>( other ).n_;
        }

    private:
        const PointCloud& fs_;
        const array::ArrayView<const double, 2> xy_;
        const array::ArrayView<const double, 1> z_;
        idx_t n_;
    };


    class IterateXYZ {
    public:
        using iterator       = IteratorXYZ;
        using const_iterator = iterator;

    public:
        IterateXYZ( const PointCloud& fs ) : fs_( fs ) {}
        iterator begin() const { return IteratorXYZ( fs_ ); }
        iterator end() const { return IteratorXYZ( fs_, false ); }

    private:
        const PointCloud& fs_;
    };

    class IteratorXY {
    public:
        IteratorXY( const PointCloud& fs, bool begin = true );

        bool next( PointXY& xyz );

        const PointXY operator*() const;

        const IteratorXY& operator++() {
            ++n_;
            return *this;
        }

        bool operator==( const IteratorXY& other ) const { return n_ == static_cast<const IteratorXY&>( other ).n_; }

        virtual bool operator!=( const IteratorXY& other ) const {
            return n_ != static_cast<const IteratorXY&>( other ).n_;
        }

    private:
        const PointCloud& fs_;
        const array::ArrayView<const double, 2> xy_;
        idx_t n_;
    };

    class IterateXY {
    public:
        using iterator       = IteratorXY;
        using const_iterator = iterator;

    public:
        IterateXY( const PointCloud& fs ) : fs_( fs ) {}
        iterator begin() const { return IteratorXY( fs_ ); }
        iterator end() const { return IteratorXY( fs_, false ); }

    private:
        const PointCloud& fs_;
    };


    class Iterate {
    public:
        Iterate( const PointCloud& fs ) : fs_( fs ) {}
        IterateXYZ xyz() const { return IterateXYZ( fs_ ); }
        IterateXY xy() const { return IterateXY( fs_ ); }

    private:
        const PointCloud& fs_;
    };

    Iterate iterate() const { return Iterate( *this ); }

private:
    Field lonlat_;
    Field vertical_;
    mutable Field ghost_;
    idx_t levels_{0};
};

//------------------------------------------------------------------------------------------------------

}  // namespace detail

//------------------------------------------------------------------------------------------------------

class PointCloud : public FunctionSpace {
public:
    PointCloud( const FunctionSpace& );
    PointCloud( const Field& points );
    PointCloud( const std::vector<PointXY>& );
    PointCloud( PointXY, const std::vector<PointXY>& );
    PointCloud( PointXYZ, const std::vector<PointXYZ>& );
    PointCloud( const Grid& grid );

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    const Field& vertical() const { return functionspace_->vertical(); }

    detail::PointCloud::Iterate iterate() const { return functionspace_->iterate(); }


private:
    const detail::PointCloud* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas
