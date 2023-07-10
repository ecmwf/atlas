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

#include <ostream>

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace helpers {

//------------------------------------------------------------------------------

class RangeBase {};

//------------------------------------------------------------------------------

class RangeFrom : public RangeBase {
public:
    RangeFrom(int start): start_(start) {}

    int start() const { return start_; }

    template <int Dim, typename View>
    int end(const View& view) const {
        return view.shape(Dim);
    }

    template <typename View>
    int end(const View& view, int i) const {
        return view.shape(i);
    }

private:
    int start_;
};

//------------------------------------------------------------------------------

class RangeTo : public RangeBase {
public:
    RangeTo(int end): end_(end) {}

    int start() const { return 0; }

    int end() const { return end_; }

private:
    int end_;
};

//------------------------------------------------------------------------------

class RangeAll : public RangeBase {
public:
    int start() const { return 0; }

    template <int Dim, typename View>
    int end(const View& view) const {
        return view.shape(Dim);
    }

    template <typename View>
    int end(const View& view, int i) const {
        return view.shape(i);
    }
    friend std::ostream& operator<<(std::ostream& out, const RangeAll& range) {
        out << ":";
        return out;
    }
};

class RangeDummy : public RangeBase {};

//------------------------------------------------------------------------------

}  // namespace helpers
#endif

//------------------------------------------------------------------------------

class Range : public helpers::RangeBase {
private:
    using From  = helpers::RangeFrom;
    using To    = helpers::RangeTo;
    using All   = helpers::RangeAll;
    using Dummy = helpers::RangeDummy;

public:
    static From from(int start) { return From(start); }
    static To to(int end) { return To(end); }
    static constexpr All all() { return All(); }
    static constexpr Dummy dummy() { return Dummy(); }

public:
    template <typename Start, typename End>
    Range(Start start, End end): start_(static_cast<int>(start)), end_(static_cast<int>(end)) {}
    int start() const { return start_; }
    int end() const { return end_; }

    Range(): Range(0, 0) {}

private:
    int start_;
    int end_;
};

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
