/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>
#include <fstream>
#include <vector>

#include "eckit/io/MemoryHandle.h"

#include "TestEnvironment.h"

namespace atlas {
namespace test {

template <typename T>
struct Matrix {
    std::vector<T> data_;
    size_t rows_;
    size_t cols_;
    size_t size() const { return rows_ * cols_; }
    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    atlas::io::DataType datatype() const { return atlas::io::make_datatype<T>(); }
    Matrix(size_t rows, size_t cols) { resize(rows, cols); }
    void resize(size_t rows, size_t cols) {
        rows_ = rows;
        cols_ = cols;
        data_.resize(size());
    }
    void assign(std::initializer_list<T> list) { data_.assign(list); }
    void assign(T value) { data_.assign(size(), value); }
};

template <typename T>
void interprete(const Matrix<T>& in, atlas::io::ArrayReference& out) {
    out = io::ArrayReference(in.data(), in.datatype(), atlas::io::ArrayShape{in.rows_, in.cols_});
}

template <typename T>
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, Matrix<T>& out) {
    atlas::io::ArrayMetadata array(metadata);
    out.resize(array.shape(0), array.shape(1));
    ::memcpy(out.data(), data, data.size());
}


struct Arrays {
    std::vector<double> v1;
    std::vector<float> v2;
    Matrix<int> v3{0, 0};
    bool operator==(const Arrays& other) const {
        return v1 == other.v1 && ::memcmp(v2.data(), other.v2.data(), v2.size() * sizeof(float)) == 0 &&
               ::memcmp(v3.data(), other.v3.data(), v3.size() * v3.datatype().size()) == 0;
    }
    bool operator!=(const Arrays& other) const { return not operator==(other); }
};

//-----------------------------------------------------------------------------

static eckit::LocalConfiguration no_compression = [] {
    eckit::LocalConfiguration c;
    c.set("compression", "none");
    return c;
}();

//-----------------------------------------------------------------------------

std::string suffix() {
    static std::string suffix = eckit::Resource<std::string>("--suffix", "");
    return suffix;
}

//-----------------------------------------------------------------------------

namespace globals {
struct TestRecord {
    Arrays data;
    TestRecord() = default;
    TestRecord(const std::function<void(Arrays&)>& initializer) { initializer(data); }
};

static TestRecord record1{[](Arrays& data) {
    data.v1 = {0, 1, 2, 3, 4};
    data.v2 = {3, 2, 1};
    data.v3.resize(3, 2);
    data.v3.assign({11, 12, 21, 22, 31, 32});
}};

static TestRecord record2{[](Arrays& data) {
    data.v1 = {0, 10, 20, 30, 40, 50};
    data.v2 = {30, 20, 10, 40};
    data.v3.resize(2, 3);
    data.v3.assign({11, 12, 13, 21, 22, 23});
}};

static TestRecord record3{[](Arrays& data) {
    data.v1.assign(1024 / 8 - 1, 2.);
    data.v2.assign(1023 * 1024 / 4 + 512 / 4, 1.);
    data.v3.resize(1024, 1024);
    data.v3.assign(3);
}};

std::vector<io::Record::URI> records;

}  // namespace globals

//-----------------------------------------------------------------------------

template <typename Length>
void write_length(Length length, const std::string& path) {
    std::ofstream file(path);
    file << length;
    file.close();
}

//-----------------------------------------------------------------------------//
//                                                                             //
//                               Writing records                               //
//                                                                             //
//-----------------------------------------------------------------------------//


CASE("Write records, each in separate file (offset=0)") {
    auto write_record = [&](const Arrays& data, const eckit::PathName& path) {
        io::RecordWriter record;
        record.set("v1", io::ref(data.v1), no_compression);
        record.set("v2", io::ref(data.v2), no_compression);
        record.set("v3", io::ref(data.v3));
        auto length = record.write(path);
        write_length(length, path + ".length");
    };

    SECTION("record1.atlas" + suffix()) {
        write_record(globals::record1.data, "record1.atlas" + suffix());
    }
    SECTION("record2.atlas" + suffix()) {
        write_record(globals::record2.data, "record2.atlas" + suffix());
    }
    SECTION("record3.atlas" + suffix()) {
        write_record(globals::record3.data, "record3.atlas" + suffix());
    }
}

//-----------------------------------------------------------------------------

CASE("Write records to same file using record.write(path,io::Mode)") {
    // This will reopen files upon every append, bad for performance

    static std::vector<size_t> lengths;
    static std::vector<size_t> offsets{0};

    auto write_record = [&](Arrays& data, const eckit::PathName& path, io::Mode mode) {
        io::RecordWriter record;
        record.set("v1", io::ref(data.v1), no_compression);
        record.set("v2", io::ref(data.v2), no_compression);
        record.set("v3", io::ref(data.v3));

        globals::records.emplace_back(io::Record::URI{path, offsets.back()});
        lengths.emplace_back(record.write(path, mode));
        offsets.emplace_back(offsets.back() + lengths.back());
    };

    SECTION("record1 -> records.atlas" + suffix()) {
        write_record(globals::record1.data, "records.atlas" + suffix(), io::Mode::write);
    }
    SECTION("record2 -> records.atlas" + suffix()) {
        write_record(globals::record2.data, "records.atlas" + suffix(), io::Mode::append);
    }
    SECTION("record3 -> records.atlas" + suffix()) {
        write_record(globals::record3.data, "records.atlas" + suffix(), io::Mode::append);
    }
}

//-----------------------------------------------------------------------------

CASE("Write records to same file using record.write(Stream) keeping Stream open") {
    // This should give exactly the same output file as previous, except no

    auto write_record = [&](const Arrays& data, io::Stream stream) {
        io::RecordWriter record;
        record.set("v1", io::ref(data.v1), no_compression);
        record.set("v2", io::ref(data.v2), no_compression);
        record.set("v3", io::ref(data.v3));

        record.write(stream);
    };

    static io::OutputFileStream stream("records.atlas" + suffix() + ".duplicate");

    SECTION("record1 -> records.atlas" + suffix() + ".duplicate") {
        EXPECT_EQ(stream.position(), globals::records[0].offset);
        write_record(globals::record1.data, stream);
    }
    SECTION("record2 -> records.atlas" + suffix() + ".duplicate") {
        EXPECT_EQ(stream.position(), globals::records[1].offset);
        write_record(globals::record2.data, stream);
    }
    SECTION("record3 -> records.atlas" + suffix() + ".duplicate") {
        EXPECT_EQ(stream.position(), globals::records[2].offset);
        write_record(globals::record3.data, stream);
    }
    SECTION("close stream") {
        stream.close();  // required because stream is a static variable
    }
}

//-----------------------------------------------------------------------------

CASE("Write master record referencing record1 and record2 and record3") {
    io::RecordWriter record;
    record.set("v1", io::link("file:record1.atlas" + suffix() + "?key=v1"));
    record.set("v2", io::link("file:record1.atlas" + suffix() + "?key=v2"));
    record.set("v3", io::link("file:record1.atlas" + suffix() + "?key=v3"));
    record.set("v4", io::link("file:record2.atlas" + suffix() + "?key=v1"));
    record.set("v5", io::link("file:record2.atlas" + suffix() + "?key=v2"));
    record.set("v6", io::link("file:record2.atlas" + suffix() + "?key=v3"));
    record.set("v7", io::link("file:record3.atlas" + suffix() + "?key=v1"));
    record.set("v8", io::link("file:record3.atlas" + suffix() + "?key=v2"));
    record.set("v9", io::link("file:record3.atlas" + suffix() + "?key=v3"));
    record.write("record.atlas" + suffix());
}

//-----------------------------------------------------------------------------

CASE("Write records in nested subdirectories") {
    auto reference_path = eckit::PathName{"atlas_test_io_refpath"};
    {
        eckit::PathName{reference_path / "links" / "1"}.mkdir();

        io::RecordWriter record;
        record.set("v1", io::ref(globals::record1.data.v1));
        record.set("v2", io::ref(globals::record1.data.v2));
        record.set("v3", io::ref(globals::record1.data.v3));
        record.set("s1", std::string("short string"));
        record.set("s2", double(1. / 3.));
        record.write(reference_path / "links" / "1" / "record.atlas" + suffix());
    }
    {
        eckit::PathName{reference_path / "links" / "2"}.mkdir();

        io::RecordWriter record;
        record.set("v1", io::ref(globals::record2.data.v1));
        record.set("v2", io::ref(globals::record2.data.v2));
        record.set("v3", io::ref(globals::record2.data.v3));
        record.set("s1", size_t(10000000000));
        record.write(reference_path / "links" / "2" / "record.atlas" + suffix());
    }
    {
        io::RecordWriter record;
        record.set("l1", io::link("file:1/record.atlas" + suffix() + "?key=v1"));
        record.set("l2", io::link("file:1/record.atlas" + suffix() + "?key=v2"));
        record.set("l3", io::link("file:1/record.atlas" + suffix() + "?key=v3"));
        record.set("l4", io::link("file:2/record.atlas" + suffix() + "?key=v1"));
        record.set("l5", io::link("file:2/record.atlas" + suffix() + "?key=v2"));
        record.set("l6", io::link("file:2/record.atlas" + suffix() + "?key=v3"));
        record.set("l7", io::link("file:1/record.atlas" + suffix() + "?key=s1"));
        record.set("l8", io::link("file:1/record.atlas" + suffix() + "?key=s2"));
        record.set("l9", io::link("file:2/record.atlas" + suffix() + "?key=s1"));
        record.write(reference_path / "links" / "record.atlas" + suffix());
    }
    {
        io::RecordWriter record;
        record.set("l1", io::link("file:links/record.atlas" + suffix() + "?key=l1"));
        record.set("l2", io::link("file:links/record.atlas" + suffix() + "?key=l2"));
        record.set("l3", io::link("file:links/record.atlas" + suffix() + "?key=l3"));
        record.set("l4", io::link("file:links/record.atlas" + suffix() + "?key=l4"));
        record.set("l5", io::link("file:links/record.atlas" + suffix() + "?key=l5"));
        record.set("l6", io::link("file:links/record.atlas" + suffix() + "?key=l6"));
        record.set("l7", io::link("file:links/record.atlas" + suffix() + "?key=l7"));
        record.set("l8", io::link("file:links/record.atlas" + suffix() + "?key=l8"));
        record.set("l9", io::link("file:links/record.atlas" + suffix() + "?key=l9"));
        record.write(reference_path / "record.atlas" + suffix());
    }
}

//-----------------------------------------------------------------------------//
//                                                                             //
//                               Reading tests                                 //
//                                                                             //
//-----------------------------------------------------------------------------//

CASE("Test RecordItemReader") {
    SECTION("file:record1.atlas" + suffix() + "?key=v2") {
        io::RecordItemReader reader{"file:record1.atlas" + suffix() + "?key=v2"};
        {
            // When we only want to read metadata
            io::Metadata metadata;
            reader.read(metadata);
            EXPECT(metadata.link() == false);
            EXPECT_EQ(metadata.type(), "array");
            EXPECT(metadata.data.compressed() == false);
            EXPECT_EQ(metadata.data.compression(), "none");
            EXPECT_EQ(metadata.data.size(), globals::record1.data.v2.size() * sizeof(float));
        }
        {
            // When we want to read both metadata and data
            io::Metadata metadata;
            io::Data data;
            reader.read(metadata, data);
            EXPECT(metadata.data.compressed() == false);
            EXPECT_EQ(metadata.data.compression(), "none");
            EXPECT_EQ(metadata.data.size(), globals::record1.data.v2.size() * sizeof(float));
            EXPECT_EQ(data.size(), metadata.data.size());
            EXPECT_EQ(data.size(), metadata.data.compressed_size());
            EXPECT(::memcmp(data, globals::record1.data.v2.data(), data.size()) == 0);
        }
    }

    SECTION("file:record.atlas" + suffix() + "?key=v9") {
        io::RecordItemReader reader{"file:record.atlas" + suffix() + "?key=v9"};
        {
            // When we only want to read metadata
            io::Metadata metadata;
            reader.read(metadata);
            EXPECT(metadata.link() == true);
            EXPECT_EQ(metadata.link().str(), "file:record3.atlas" + suffix() + "?key=v3");
            EXPECT_EQ(metadata.type(), "array");
        }
        {
            // When we want to read both metadata and data
            io::Metadata metadata;
            io::Data data;
            reader.read(metadata, data);
            EXPECT(metadata.data.compressed() == (io::defaults::compression_algorithm() != "none"));
            EXPECT_EQ(metadata.data.compression(), io::defaults::compression_algorithm());
            EXPECT_EQ(metadata.data.size(), globals::record3.data.v3.size() * sizeof(int));
            EXPECT_EQ(data.size(), metadata.data.compressed_size());
        }
    }
}

//-----------------------------------------------------------------------------

CASE("Read records from different files") {
    Arrays data1, data2, data3;

    auto read_record = [&](const eckit::PathName& path, Arrays& data) {
        io::RecordReader record(path);
        record.read("v1", data.v1).wait();
        record.read("v2", data.v2).wait();
        record.read("v3", data.v3).wait();
    };

    read_record("record1.atlas" + suffix(), data1);
    read_record("record2.atlas" + suffix(), data2);
    read_record("record3.atlas" + suffix(), data3);

    EXPECT(data1 == globals::record1.data);
    EXPECT(data2 == globals::record2.data);
    EXPECT(data3 == globals::record3.data);
}

//-----------------------------------------------------------------------------

CASE("Read multiple records from same file") {
    Arrays data1, data2;
    io::RecordReader record1(globals::records[0]);
    io::RecordReader record2(globals::records[1]);

    record1.read("v1", data1.v1).wait();
    record1.read("v2", data1.v2).wait();
    record1.read("v3", data1.v3).wait();

    record2.read("v1", data2.v1).wait();
    record2.read("v2", data2.v2).wait();
    record2.read("v3", data2.v3).wait();

    EXPECT(data1 == globals::record1.data);
    EXPECT(data2 == globals::record2.data);
}

//-----------------------------------------------------------------------------

CASE("Write master record referencing record1 and record2") {
    io::RecordWriter record;
    record.set("v1", io::link("file:record1.atlas" + suffix() + "?key=v1"));
    record.set("v2", io::link("file:record1.atlas" + suffix() + "?key=v2"));
    record.set("v3", io::link("file:record1.atlas" + suffix() + "?key=v3"));
    record.set("v4", io::link("file:record2.atlas" + suffix() + "?key=v1"));
    record.set("v5", io::link("file:record2.atlas" + suffix() + "?key=v2"));
    record.set("v6", io::link("file:record2.atlas" + suffix() + "?key=v3"));
    record.write("record.atlas" + suffix());
}


//-----------------------------------------------------------------------------

CASE("Read master record") {
    Arrays data1, data2;
    io::RecordReader record("record.atlas" + suffix());

    eckit::Log::info() << "record.metadata(\"v1\"): " << record.metadata("v1") << std::endl;


    record.read("v1", data1.v1).wait();
    record.read("v2", data1.v2).wait();
    record.read("v3", data1.v3).wait();
    record.read("v4", data2.v1).wait();
    record.read("v5", data2.v2).wait();
    record.read("v6", data2.v3).wait();

    EXPECT(data1 == globals::record1.data);
    EXPECT(data2 == globals::record2.data);
}

//-----------------------------------------------------------------------------

CASE("Async read") {
    Arrays data1, data2;
    io::RecordReader record("record.atlas" + suffix());

    // Request reads
    record.read("v1", data1.v1);
    record.read("v2", data1.v2);
    record.read("v3", data1.v3);
    record.read("v4", data2.v1);
    record.read("v5", data2.v2);
    record.read("v6", data2.v3);

    // Wait for specific requests
    record.wait("v4");
    record.wait("v5");
    record.wait("v6");

    // Should have completed
    EXPECT(data2 == globals::record2.data);

    // Should not be complete yet
    EXPECT(data1 != globals::record1.data);

    // Wait for all requests;
    record.wait();

    // Should have completed
    EXPECT(data1 == globals::record1.data);
}

//-----------------------------------------------------------------------------

CASE("Recursive Write/read records in nested subdirectories") {
    auto reference_path = eckit::PathName{"atlas_test_io_refpath"};

    //  Read

    Arrays data1, data2;
    io::RecordReader record(reference_path / "record.atlas" + suffix());

    record.read("l1", data1.v1).wait();
    record.read("l2", data1.v2).wait();
    record.read("l3", data1.v3).wait();

    record.read("l4", data2.v1).wait();
    record.read("l5", data2.v2).wait();
    record.read("l6", data2.v3).wait();

    std::string l7;
    double l8;
    size_t l9;
    record.read("l7", l7).wait();
    record.read("l8", l8).wait();
    record.read("l9", l9).wait();

    EXPECT(data1 == globals::record1.data);
    EXPECT(data2 == globals::record2.data);
    EXPECT_EQ(l7, "short string");
    EXPECT_EQ(l8, 1. / 3.);
    EXPECT_EQ(l9, 10000000000ul);
}

//-----------------------------------------------------------------------------

CASE("Write record to memory") {
    const auto& data_write = globals::record3.data;
    const auto& v1         = globals::record3.data.v1;
    const auto& v2         = globals::record3.data.v2;
    const auto& v3         = globals::record3.data.v3;

    eckit::Buffer memory;

    // write
    {
        ATLAS_IO_TRACE("write");
        io::RecordWriter record;
        record.compression(false);
        record.checksum(false);
        record.set("v1", io::ref(v1));
        record.set("v2", io::ref(v2));
        record.set("v3", io::ref(v3));

        memory.resize(record.estimateMaximumSize());

        eckit::Log::info() << "memory.size() : " << memory.size() << std::endl;
        ;

        eckit::MemoryHandle datahandle_out{memory};
        datahandle_out.openForWrite(0);
        auto record_length = record.write(datahandle_out);
        datahandle_out.close();

        // Without compression, this should be exact
        EXPECT_EQ(memory.size(), record_length);
    }

    // read with individual RecordItemReader
    {
        ATLAS_IO_TRACE("read with RecordItemReader");

        io::Session session;

        eckit::MemoryHandle datahandle_in{memory};
        datahandle_in.openForRead();

        {
            io::RecordItemReader reader(datahandle_in, "v1");
            io::Metadata metadata;
            io::Data data;
            reader.read(metadata, data);
            EXPECT(::memcmp(data, data_write.v1.data(), data.size()) == 0);
        }
        {
            io::RecordItemReader reader(datahandle_in, "v2");
            io::Metadata metadata;
            io::Data data;
            reader.read(metadata, data);
            EXPECT(::memcmp(data, data_write.v2.data(), data.size()) == 0);
        }
        {
            io::RecordItemReader reader(datahandle_in, "v3");
            io::Metadata metadata;
            io::Data data;
            reader.read(metadata, data);
            EXPECT(::memcmp(data, data_write.v3.data(), data.size()) == 0);
        }
        datahandle_in.close();
    }

    // read with RecordReader
    {
        ATLAS_IO_TRACE("read with RecordReader");
        Arrays data_read;

        eckit::MemoryHandle datahandle_in{memory};
        datahandle_in.openForRead();

        io::RecordReader reader(datahandle_in);
        reader.read("v1", data_read.v1);
        reader.read("v2", data_read.v2);
        reader.read("v3", data_read.v3);
        reader.wait();

        datahandle_in.close();

        EXPECT(data_read == data_write);
    }
}

//-----------------------------------------------------------------------------//
//                                                                             //
//                               Reading tests                                 //
//                                                                             //
//-----------------------------------------------------------------------------//

CASE("RecordPrinter") {
    SECTION("table") {
        eckit::LocalConfiguration table_with_details;
        table_with_details.set("format", "table");
        table_with_details.set("details", true);

        io::RecordPrinter record{eckit::PathName("record1.atlas" + suffix()), table_with_details};
        std::stringstream out;
        EXPECT_NO_THROW(out << record);
        eckit::Log::debug() << out.str();
    }

    SECTION("yaml") {
        eckit::LocalConfiguration yaml_with_details;
        yaml_with_details.set("format", "yaml");
        yaml_with_details.set("details", true);

        io::RecordPrinter record{eckit::PathName("record1.atlas" + suffix()), yaml_with_details};
        std::stringstream out;
        EXPECT_NO_THROW(out << record);
        eckit::Log::debug() << out.str();
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
