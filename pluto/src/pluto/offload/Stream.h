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

#include <memory>
#include <functional>


namespace pluto {

class Stream {
    using stream_t = void*;
public:
    explicit Stream(stream_t stream);
    Stream();

    [[nodiscard]] stream_t value() const {
        return *stream_;
    }

    template <typename deviceStream_t>
    [[nodiscard]] deviceStream_t value() const {
        return reinterpret_cast<deviceStream_t>(*stream_);
    }

    void wait() const;
private:
    std::unique_ptr<stream_t, std::function<void(stream_t*)>> stream_;
};


const Stream& default_stream();

const Stream& get_default_stream();

void set_default_stream(const Stream&);

void wait(const Stream&);

class [[nodiscard]] scoped_default_stream {
public:

    scoped_default_stream(const Stream& s) :
        saved_(get_default_stream()) {
        set_default_stream(s);
    }

    ~scoped_default_stream() {
        set_default_stream(saved_);
    }
private:
   const Stream& saved_;
};

}
