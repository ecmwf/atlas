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

class stream {
    using stream_t = void*;
public:
    explicit stream(stream_t);
    stream();

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


const stream& default_stream();

const stream& get_current_stream();

void set_current_stream(const stream&);

void wait(const stream&);

class [[nodiscard]] current_stream {
public:

    current_stream(const stream& s) :
        saved_(get_current_stream()) {
        set_current_stream(s);
    }

    ~current_stream() {
        set_current_stream(saved_);
    }
private:
   const stream& saved_;
};

}
