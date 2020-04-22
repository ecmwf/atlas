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
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"


namespace atlas {
class Grid;
}


namespace atlas {
namespace functionspace {


namespace detail {

class Points final : public FunctionSpaceImpl {
public:
    Points( const Grid& );
    ~Points() override;

    using FunctionSpaceImpl::createField;
    Field createField( const Field&, const eckit::Configuration& ) const override;
    Field createField( const eckit::Configuration& ) const override;

    void haloExchange( const FieldSet&, bool /*on_device*/ = false ) const override;
    void haloExchange( const Field&, bool /*on_device*/ = false ) const override;

    idx_t size() const override;
    size_t footprint() const override;
    std::string distribution() const override;
    std::string type() const override;

    atlas::Field lonlat() const override { return lonlat_; }
    atlas::Field ghost() const override;


    template <typename T>
    struct IteratorT {
        IteratorT( const Field& field, bool begin ) :
            lonlat_( field ),
            view_( array::make_view<double, 2>( lonlat_ ) ),
            m_( lonlat_.shape( 0 ) ),
            n_( begin ? 0 : m_ ) {}

        bool next( T& p ) {
            if ( n_ < m_ ) {
                p = operator*();
                ++n_;
                return true;
            }
            return false;
        }

        const IteratorT& operator++() {
            ++n_;
            return *this;
        }

        const T operator*() const { ATLAS_NOTIMPLEMENTED; }

        bool operator==( const IteratorT& other ) const { return n_ == other.n_; }
        bool operator!=( const IteratorT& other ) const { return n_ != other.n_; }

    private:
        const Field lonlat_;
        const array::ArrayView<const double, 2> view_;
        idx_t m_;
        idx_t n_;
    };


    template <typename T>
    struct IterateT {
        using iterator       = IteratorT<T>;
        using const_iterator = iterator;

        IterateT( const Field& lonlat ) : lonlat_( lonlat ) {}
        iterator begin() const { return {lonlat_, true}; }
        iterator end() const { return {lonlat_, false}; }

    private:
        const Field lonlat_;
    };


    struct Iterate {
        Iterate( const Points& points ) : lonlat_( points.lonlat() ) {}
        IterateT<PointXYZ> xyz() const { return {lonlat_}; }
        IterateT<PointLonLat> lonlat() const { return {lonlat_}; }

    private:
        const Field lonlat_;
    };

    Iterate iterate() const { return Iterate( *this ); }

private:
    Field lonlat_;
    mutable Field ghost_;
};


}  // namespace detail


class Points : public FunctionSpace {
public:
    Points( const Grid& );
    Points( const FunctionSpace& );

    operator bool() const { return functionspace_ != nullptr; }

    Field lonlat() const { return functionspace_->lonlat(); }
    Field ghost() const { return functionspace_->ghost(); }

    detail::Points::Iterate iterate() const { return functionspace_->iterate(); }

private:
    const detail::Points* functionspace_;
};


}  // namespace functionspace
}  // namespace atlas
