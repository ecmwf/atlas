/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "stream.h"

#include "hic/hic.h"
#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/trace.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

static void* stream0_underlying_ = nullptr;

#if PLUTO_HAVE_HIC
stream::stream(stream_t& s):
    stream_{&s, [](stream_t*) {
                // wrapping, no delete
            }} {}

stream::stream():
    stream_{[]() {
                if (devices()) {
                    hicStream_t* s = new hicStream_t;
                    HIC_CALL(hicStreamCreate(s));
                    return reinterpret_cast<stream_t*>(s);
                }
                return &stream0_underlying_;
            }(),
            [](stream_t* s) {
                if (s != &stream0_underlying_) {
                    HIC_CALL(hicStreamDestroy(*reinterpret_cast<hicStream_t*>(s)));
                    delete s;
                }
            }} {}

void stream::wait() const {
    if constexpr (LOG) {
        std::cout << "               = hicStreamSynchronize(stream:" << value() << ")" << std::endl;
    }
    HIC_CALL(hicStreamSynchronize(value<hicStream_t>()));
}

void stream_view::wait() const {
    if constexpr (LOG) {
        std::cout << "               = hicStreamSynchronize(stream:" << value() << ")" << std::endl;
    }
    HIC_CALL(hicStreamSynchronize(value<hicStream_t>()));
}


#else
stream::stream(stream_t& stream):
    stream_{&stream, [](stream_t*) {
                // wrapping, no delete
            }} {}

stream::stream(): stream(stream0_underlying_) {}

void stream::wait() const {
    // Nothing
}

void stream_view::wait() const {
    // Nothing
}

#endif

static stream stream0_{stream0_underlying_};

stream_view default_stream() {
    return stream_view{stream0_};
}

static thread_local stream_view current_stream_ = default_stream();


stream_view get_stream() {
    return current_stream_;
}

void set_stream(stream_view s) {
    current_stream_ = s;
}

void wait(stream_view s) {
    s.wait();
}


}  // namespace pluto
