/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Stream.h"

#include "pluto/pluto_config.h"
#include "hic/hic.h"
#include "pluto/util/Trace.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

#if PLUTO_HAVE_HIC
Stream::Stream(stream_t stream)
    : stream_ {
        &stream,
        [](stream_t* stream) {
            // wrapping, no delete
        }
    } {
}

Stream::Stream()
    :stream_{
        []() {
            hicStream_t* stream = new hicStream_t;
            hicStreamCreate(stream);
            return reinterpret_cast<stream_t*>(stream);
        }(),
        [](stream_t* stream) {
            hicStreamDestroy(*reinterpret_cast<hicStream_t*>(stream));
            delete stream;
        }
    } {
}

void Stream::wait() const {
    if constexpr(LOG) {
        std::cout << "               = hicStreamSynchronize(stream:"<<value()<<")" << std::endl;
    }
    hicStreamSynchronize(value<hicStream_t>());
}

#else
Stream::Stream(stream_t stream)
    : stream_ {
        &stream,
        [](stream_t* stream) {
            // wrapping, no delete
        }
    } {
}

Stream::Stream()
    :stream_{
        []() {
            static int stream = 0;
            return reinterpret_cast<stream_t*>(&stream);
        }(),
        [](stream_t* stream) {
            // No deletion of wrapped static variable
        }
    } {
}

void Stream::wait() const {
    // Nothing
}
#endif

static Stream stream0_{nullptr};

static const Stream* default_stream_ = &stream0_;

const Stream& get_default_stream() {
    return *default_stream_;
}

void set_default_stream(const Stream& s) {
    default_stream_ = &s;
}

void wait(const Stream& stream) {
    stream.wait();
}

}
