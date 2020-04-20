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

    const atlas::Field xyz() const { return xyz_; }


    template <typename T>
    struct IteratorT {
        IteratorT( const Field& field, bool begin = true ) :
            field_( field ),
            view_( array::make_view<double, 2>( field_ ) ),
            m_( field_.shape( 0 ) ),
            n_( begin ? 0 : m_ ) {}

        bool next( PointXYZ& p ) {
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
        const Field& field_;
        const array::ArrayView<const double, 2> view_;
        idx_t m_;
        idx_t n_;
    };


    template <typename T>
    struct IterateT {
        using iterator       = IteratorT<T>;
        using const_iterator = iterator;

        IterateT( const Field& field ) : field_( field ) {}
        iterator begin() const { return {field_}; }
        iterator end() const { return {field_, false}; }

    private:
        const Field& field_;
    };


    struct Iterate {
        Iterate( const Points& fs ) : fs_( fs ) {}
        IterateT<PointXYZ> xyz() const { return {fs_.xyz()}; }
        IterateT<PointLonLat> lonlat() const { return {fs_.lonlat()}; }

    private:
        const Points& fs_;
    };

    Iterate iterate() const { return Iterate( *this ); }

private:
    Field lonlat_;
    Field xyz_;
    mutable Field ghost_;
};


}  // namespace detail


class Points : public FunctionSpace {
public:
    Points( const Grid& );
    Points( const FunctionSpace& );

    operator bool() const { return functionspace_ != nullptr; }

    Field lonlat() const { return functionspace_->lonlat(); }
    Field xyz() const { return functionspace_->xyz(); }
    Field ghost() const { return functionspace_->ghost(); }

    detail::Points::Iterate iterate() const { return functionspace_->iterate(); }

private:
    const detail::Points* functionspace_;
};


}  // namespace functionspace
}  // namespace atlas
