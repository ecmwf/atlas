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
#include "pluto/trace.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

#if PLUTO_HAVE_HIC
stream::stream(stream_t s):
    stream_{&s, [](stream_t*) {
                // wrapping, no delete
            }} {}

stream::stream():
    stream_{[]() {
                hicStream_t* s = new hicStream_t;
                HIC_CALL(hicStreamCreate(s));
                return reinterpret_cast<stream_t*>(s);
            }(),
            [](stream_t* s) {
                HIC_CALL(hicStreamDestroy(*reinterpret_cast<hicStream_t*>(s)));
                delete s;
            }} {}

void stream::wait() const {
    if constexpr (LOG) {
        std::cout << "               = hicStreamSynchronize(stream:" << value() << ")" << std::endl;
    }
    HIC_CALL(hicStreamSynchronize(value<hicStream_t>()));
}

#else
stream::stream(stream_t stream):
    stream_{&stream, [](stream_t* stream) {
                // wrapping, no delete
            }} {}

stream::stream():
    stream_{[]() {
                static int s = 0;
                return reinterpret_cast<stream_t*>(&s);
            }(),
            [](stream_t*) {
                // No deletion of wrapped static variable
            }} {}

void stream::wait() const {
    // Nothing
}
#endif

static stream stream0_{nullptr};

static const stream* default_stream_ = &stream0_;

const stream& default_stream() {
    return stream0_;
}

const stream& get_current_stream() {
    return *default_stream_;
}

void set_current_stream(const stream& s) {
    default_stream_ = &s;
}

void wait(const stream& s) {
    s.wait();
}

}  // namespace pluto
