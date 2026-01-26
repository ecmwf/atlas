/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <functional>
#include <memory>


namespace pluto {

class stream {
    using stream_t = void*;

public:
    explicit stream(stream_t&);
    stream();

    [[nodiscard]] stream_t value() const { return *stream_; }

    template <typename deviceStream_t>
    [[nodiscard]] deviceStream_t value() const {
        return reinterpret_cast<deviceStream_t>(*stream_);
    }

    void wait() const;

    const stream_t* data() const { return stream_.get(); }

private:
    std::unique_ptr<stream_t, std::function<void(stream_t*)>> stream_;
};

class stream_view {
    using stream_t = void*;

public:
    stream_view(const stream& s): stream_(s.data()) {}
    stream_view() {}

    bool empty() const { return stream_ == nullptr; }

    [[nodiscard]] const stream_t* data() const { return stream_; }

    [[nodiscard]] stream_t value() const { return *stream_; }

    template <typename deviceStream_t>
    [[nodiscard]] const deviceStream_t& value() const {
        return reinterpret_cast<const deviceStream_t&>(*stream_);
    }

    void wait() const;

private:
    const stream_t* stream_{nullptr};
};


stream_view default_stream();

stream_view get_stream();

void set_stream(stream_view);

void wait(stream_view);

class [[nodiscard]] scoped_stream {
public:
    scoped_stream(stream_view s): saved_(get_stream()) { set_stream(s); }

    ~scoped_stream() { set_stream(saved_); }

private:
    stream_view saved_;
};


}  // namespace pluto
