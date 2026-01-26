/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <sstream>

#include "pluto/pluto.h"
#include "atlas/field.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

CASE("default resource") {
    std::stringstream pluto_trace_stream;
    
    pluto::trace::enable(true);
    pluto::trace::set(pluto_trace_stream);
    Field field_1("field_1", make_datatype<double>(), array::make_shape(20, 5));
    {
        Field field_2("field_2", make_datatype<double>(), array::make_shape(200, 5));
    }
    Field field_3("field_3", make_datatype<double>(), array::make_shape(20, 10));

    EXPECT_EQ( pluto::memory::host.total_allocations() , 3 );
    EXPECT_EQ( pluto::memory::host.allocations() , 2 );
    EXPECT_EQ( pluto::memory::host.bytes() , 2400 );
    EXPECT_EQ( pluto::memory::host.high_watermark() , 8800 );

    auto pluto_trace = pluto_trace_stream.str();
    std::cout << "pluto_trace:\n" << pluto_trace << std::endl;

    auto pluto_trace_contains = [&pluto_trace](const std::string& str) {
        return (pluto_trace.find(str) != std::string::npos);
    };

    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::host_resource::allocate(label=field_1, bytes=800.00B, alignment=256)"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::host_resource::allocate(label=field_2, bytes=7.81K, alignment=256)"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::host_resource::deallocate(label=field_2"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::host_resource::allocate(label=field_3, bytes=1.56K, alignment=256)"));

    pluto::trace::enable(false);
}

CASE("pinned_pool_resource") {

    pluto::host::set_default_resource("pluto::pinned_pool_resource");

    std::stringstream pluto_trace_stream;
    
    pluto::trace::enable(true);
    pluto::trace::set(pluto_trace_stream);
    Field field_1("field_1", make_datatype<double>(), array::make_shape(20, 5));
    {
        Field field_2("field_2", make_datatype<double>(), array::make_shape(200, 5));
    }
    Field field_3("field_3", make_datatype<double>(), array::make_shape(20, 10));

    EXPECT_EQ( pluto::memory::pinned.total_allocations() , 1 );
    EXPECT_EQ( pluto::memory::pinned.allocations() , 1 );
    EXPECT_EQ( pluto::memory::pinned.bytes() , 256*1024*1024 ); // 256 mb
    EXPECT_EQ( pluto::memory::pinned.high_watermark() , 256*1024*1024 );

    EXPECT_EQ( pluto::memory::pinned_pool.total_allocations() , 3 );
    EXPECT_EQ( pluto::memory::pinned_pool.allocations() , 2 );
    EXPECT_EQ( pluto::memory::pinned_pool.bytes() , 2400 );
    EXPECT_EQ( pluto::memory::pinned_pool.high_watermark() , 8800 );

    auto pluto_trace = pluto_trace_stream.str();
    std::cout << "pluto_trace:\n" << pluto_trace << std::endl;

    auto pluto_trace_contains = [&pluto_trace](const std::string& str) {
        return (pluto_trace.find(str) != std::string::npos);
    };

    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::pinned_resource::allocate(label=pool_chunk, bytes=256.00M, alignment=256)"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::pinned_pool_resource::allocate(label=field_1, bytes=800.00B, alignment=256)"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::pinned_pool_resource::allocate(label=field_2, bytes=7.81K, alignment=256)"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::pinned_pool_resource::deallocate(label=field_2"));
    EXPECT(pluto_trace_contains("PLUTO_TRACE pluto::pinned_pool_resource::allocate(label=field_3, bytes=1.56K, alignment=256)"));

    pluto::trace::enable(false);

}


}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
