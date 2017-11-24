/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

namespace helpers {

//------------------------------------------------------------------------------

class RangeBase {};

//------------------------------------------------------------------------------

class RangeFrom : public RangeBase {
public:
  RangeFrom(int start) : start_(start) {}

  int start() const { return start_; }

  template < int Dim, typename View >
  int end(const View& view) const { return view.shape(Dim); }

private:
  int start_;
};

//------------------------------------------------------------------------------

class RangeTo : public RangeBase {
public:
  RangeTo(int end) : end_(end) {}

  int start() const { return 0; }

  int end() const { return end_; }

private:
  int end_;
};

//------------------------------------------------------------------------------

class RangeAll : public RangeBase {
public:

  int start() const { return 0; }

  template < int Dim, typename View >
  int end(const View& view) const { return view.shape(Dim); }
};

//------------------------------------------------------------------------------

} // helpers

//------------------------------------------------------------------------------

class Range : public helpers::RangeBase{
private:
  using From = helpers::RangeFrom;
  using To   = helpers::RangeTo;
  using All  = helpers::RangeAll;

public:

  static From from(int start) { return From(start); }
  static To   to(int end) { return To(end); }
  static All  all() { return All(); }

public:

  Range(int start, int end ) : start_(start), end_(end) {}
  int start() const { return start_; }
  int end() const { return end_; }

private:
  int start_;
  int end_;
};

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
