/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/runtime/Memory.h"
#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

struct CustomMemoryResource : public pluto::memory_resource {
    void* do_allocate(std::size_t size, std::size_t alignment) override {
        std::cout << "   + custom allocate" << std::endl;
        return pluto::new_delete_resource()->allocate(size, alignment);
    }
    void do_deallocate(void* ptr, std::size_t size, std::size_t alignment) override {
        std::cout << "   - custom deallocate" << std::endl;
        pluto::new_delete_resource()->deallocate(ptr, size, alignment);
    }
    bool do_is_equal(const memory_resource& other) const override {
        return true;
    }
};

void run_allocator(pluto::allocator<double>&& allocator) {
    std::size_t size = 10;
    std::cout << "   label = " << memory::label::get() << std::endl;
    auto resource = pluto::get_registered_name(allocator.resource());
    std::cout << "   resource = " << resource << std::endl;
    double* data = allocator.allocate(memory::label::get(),size);
    allocator.deallocate(memory::label::get(),data, size);
}

void run_default() {
    std::cout << "run_default()" << std::endl;
    memory::label label("default");
    run_allocator(pluto::allocator<double>());
}

void run_resource(pluto::memory_resource* mr) {
    std::cout << "run_resource()" << std::endl;
    memory::label label("pluto::memory_resource*");
    run_allocator(pluto::host::allocator<double>(mr));
}

void run_registered(std::string_view resource) {
    std::cout << "run_registered("<<resource<<")" << std::endl;
    memory::label label(resource);
    run_resource(pluto::get_registered_resource(resource));
}

void run_scoped_default(std::string_view resource) {
    std::cout << "run_scoped_default("<<resource<<")" << std::endl;
    atlas::memory::scope mem_scope;
    pluto::host::set_default_resource(resource);
    run_default();
}

CASE("test scope") {
    atlas::memory::set_unified(true);
    EXPECT_EQ(atlas::memory::get_unified(),true);
    atlas::memory::set_unified(false);
    EXPECT_EQ(atlas::memory::get_unified(),false);
    atlas::memory::scope::push();
    {
        EXPECT_EQ(atlas::memory::get_unified(),false);
        atlas::memory::set_unified(true);
        EXPECT_EQ(atlas::memory::get_unified(),true);
    }
    atlas::memory::scope::pop();
    EXPECT_EQ(atlas::memory::get_unified(),false);

    // Now nested scope
    atlas::memory::scope::push();
    {
        atlas::memory::set_unified(true);
        atlas::memory::scope::push();
        {
            EXPECT_EQ(atlas::memory::get_unified(),true);
            atlas::memory::set_unified(false);
            atlas::memory::scope::push();
            {
                EXPECT_EQ(atlas::memory::get_unified(),false);
                atlas::memory::set_unified(true);
                EXPECT_EQ(atlas::memory::get_unified(),true);
            }
            atlas::memory::scope::pop();
            EXPECT_EQ(atlas::memory::get_unified(),false);
        }
        atlas::memory::scope::pop();
        EXPECT_EQ(atlas::memory::get_unified(),true);
    }
    atlas::memory::scope::pop();
    EXPECT_EQ(atlas::memory::get_unified(),false); 
}

CASE("test scope alloc") {
    run_default();
    run_scoped_default("pluto::pinned_resource");
    run_scoped_default("pluto::new_delete_resource");
    run_scoped_default("pluto::pinned_pool_resource");
    pluto::pinned_pool_resource()->release();
}

// --------------------------------------------------------------------------


CASE("test extension") {
    CustomMemoryResource mr;
    pluto::register_resource("custom_resource", &mr);

    run_scoped_default("custom_resource");

    pluto::unregister_resource("custom_resource");
}

// --------------------------------------------------------------------------

CASE("test context") {
    CustomMemoryResource mr;
    pluto::Register mr_register("custom_resource", &mr);

    memory::scope::push();
    pluto::host::set_default_resource("custom_resource");
    memory::set_unified(false);
    memory::register_context("custom");
    memory::scope::pop();

    run_default();

    memory::scope::push();
    pluto::host::set_default_resource("pluto::managed_resource");
    memory::set_unified(true);
    memory::register_context("unified");
    memory::scope::pop();

    EXPECT(memory::context_exists("custom"));
    EXPECT(memory::context_exists("unified"));

    run_default();

    memory::scope::push();
    memory::set_context("custom");
    run_default();

    memory::set_context("unified");
    run_default();
    memory::scope::pop();

    memory::unregister_context("custom");
    memory::unregister_context("unified");

}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
