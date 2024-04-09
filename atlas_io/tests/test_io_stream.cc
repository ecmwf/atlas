/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_io/atlas-io.h"

#include "eckit/io/FileHandle.h"
#include "eckit/io/PooledHandle.h"
#include "TestEnvironment.h"


namespace atlas {
namespace test {

CASE("Stream interoperability with eckit::DataHandle") {
    SECTION("own pointer") {
        io::Stream s;
        {
            eckit::DataHandle* datahandle = new eckit::FileHandle("test_io_session.data");
            datahandle->openForWrite(0);
            s = io::Stream{datahandle};
        }
        s.datahandle().close();
    }
    SECTION("shared pointer") {
        io::Stream s;
        {
            std::shared_ptr<eckit::DataHandle> datahandle = std::make_shared<eckit::FileHandle>("test_io_session.data");
            datahandle->openForWrite(0);
            s = io::Stream{datahandle};
        }
        s.datahandle().close();
    }
    SECTION("reference") {
        io::Stream s;
        eckit::FileHandle datahandle("test_io_session.data");
        datahandle.openForWrite(0);
        s = io::Stream{datahandle};
        s.datahandle().close();
    }
}

CASE("Test seek-for-write works when opening OutputFileStream for append") {
    std::string s1("write    \n");
    std::string s2("append   \n");
    std::string s3("overwrite\n");
    {
        ATLAS_IO_TRACE("write");
        io::Stream f = io::OutputFileStream("append-test");
        f.write(s1.c_str(), s1.size());
    }
    {
        ATLAS_IO_TRACE("append");
        io::Stream f = io::OutputFileStream("append-test", io::Mode::append);
        auto offset  = f.position();
        f.write(s2.c_str(), s2.size());

        // Rewind to beginning of append
        f.seek(offset);
        f.write(s3.c_str(), s3.size());
    }
    {
        ATLAS_IO_TRACE("read");
        io::Stream f         = io::InputFileStream("append-test");
        std::string expected = s1 + s3;
        std::string read(expected.size(), ' ');
        f.read(const_cast<char*>(read.data()), read.size());
        EXPECT_EQ(read, expected);
    }
}


CASE("Opening same file in same scope") {
    // Opening same file within same scope will avoid opening it multiple times, good for perfmance

    // write a file
    {
        io::OutputFileStream out("test_io_session.data");
        out.write("line1\n", 6);
        out.write("line2\n", 6);
        out.write("line3\n", 6);
    }

    std::string l1(5, ' '), l2(5, ' '), l3(5, ' ');

    io::Stream f1 = io::InputFileStream{"test_io_session.data"};
    f1.seek(0 * 6);
    f1.read(const_cast<char*>(l1.data()), 5);

    io::Stream f2 = io::InputFileStream{"test_io_session.data"};
    f2.seek(1 * 6);
    f2.read(const_cast<char*>(l2.data()), 5);

    io::Stream f3 = io::InputFileStream{"test_io_session.data"};
    f3.seek(2 * 6);
    f3.read(const_cast<char*>(l3.data()), 5);

    EXPECT_EQ(l1, "line1");
    EXPECT_EQ(l2, "line2");
    EXPECT_EQ(l3, "line3");

    auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f1.datahandle());
    EXPECT_EQ(pooled_handle.nbOpens(), 1);
    EXPECT_EQ(pooled_handle.nbSeeks(), 3);
    EXPECT_EQ(pooled_handle.nbReads(), 3);
}

CASE("Opening same file in parallel scopes") {
    // Files are opened and closed within each scope, bad for performance

    // write a file
    {
        io::OutputFileStream out("test_io_session.data");
        out.write("line1\n", 6);
        out.write("line2\n", 6);
        out.write("line3\n", 6);
    }

    std::string l1(5, ' '), l2(5, ' '), l3(5, ' ');

    {
        io::Stream f1 = io::InputFileStream{"test_io_session.data"};
        f1.seek(0 * 6);
        f1.read(const_cast<char*>(l1.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f1.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 1);
        EXPECT_EQ(pooled_handle.nbReads(), 1);
    }
    {
        io::Stream f2 = io::InputFileStream{"test_io_session.data"};
        f2.seek(1 * 6);
        f2.read(const_cast<char*>(l2.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f2.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 1);
        EXPECT_EQ(pooled_handle.nbReads(), 1);
    }
    {
        io::Stream f3 = io::InputFileStream{"test_io_session.data"};
        f3.seek(2 * 6);
        f3.read(const_cast<char*>(l3.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f3.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 1);
        EXPECT_EQ(pooled_handle.nbReads(), 1);
    }
}

CASE("Opening same file in parallel scopes with Session") {
    // Declaring this in an outer scope will keep storage of InputFileStream
    // within nested scopes, so that files will not be opened/closed repeatedly

    io::Session session;

    // write a file
    {
        io::OutputFileStream out("test_io_session.data");
        out.write("line1\n", 6);
        out.write("line2\n", 6);
        out.write("line3\n", 6);
    }


    std::string l1(5, ' '), l2(5, ' '), l3(5, ' ');

    {
        io::Stream f1 = io::InputFileStream{"test_io_session.data"};
        f1.seek(0 * 6);
        f1.read(const_cast<char*>(l1.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f1.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 1);
        EXPECT_EQ(pooled_handle.nbReads(), 1);
    }
    {
        io::Stream f2 = io::InputFileStream{"test_io_session.data"};
        f2.seek(1 * 6);
        f2.read(const_cast<char*>(l2.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f2.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 2);
        EXPECT_EQ(pooled_handle.nbReads(), 2);
    }
    {
        io::Stream f3 = io::InputFileStream{"test_io_session.data"};
        f3.seek(2 * 6);
        f3.read(const_cast<char*>(l3.data()), 5);
        auto& pooled_handle = dynamic_cast<eckit::PooledHandle&>(f3.datahandle());
        EXPECT_EQ(pooled_handle.nbOpens(), 1);
        EXPECT_EQ(pooled_handle.nbSeeks(), 3);
        EXPECT_EQ(pooled_handle.nbReads(), 3);
    }
}


}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
