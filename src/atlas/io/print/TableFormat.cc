/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "TableFormat.h"

#include <iomanip>

#include "atlas/io/Exceptions.h"
#include "atlas/io/Metadata.h"
#include "atlas/io/Record.h"
#include "atlas/io/RecordItemReader.h"
#include "atlas/io/Session.h"
#include "atlas/io/print/Bytes.h"
#include "atlas/io/types/array/ArrayReference.h"
#include "atlas/io/types/scalar.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace io {


std::ostream& operator<<( std::ostream& out, const MetadataPrettyPrintBase& p ) {
    p.print( out );
    return out;
}

std::string MetadataPrettyPrintBase::str() const {
    std::stringstream s;
    print( s );
    return s.str();
}

class DefaultMetadataPrettyPrint : public MetadataPrettyPrintBase {
public:
    DefaultMetadataPrettyPrint( const Metadata& ) {}
    void print( std::ostream& ) const override {}
};

class ArrayMetadataPrettyPrint : public MetadataPrettyPrintBase {
public:
    ArrayMetadataPrettyPrint( const Metadata& m ) : metadata_( m ) {}
    void print( std::ostream& out ) const override {
        std::string type = metadata_.getString( "type" );
        ATLAS_ASSERT( type == "array" );
        ArrayMetadata array( metadata_ );
        out << std::setw( 7 ) << std::left << array.datatype().str();
        out << "[";
        for ( int i = 0; i < array.rank(); ++i ) {
            out << array.shape( i );
            if ( i < array.rank() - 1 ) {
                out << ",";
            }
        }
        out << "]";
    }

private:
    Metadata metadata_;
};

class StringMetadataPrettyPrint : public MetadataPrettyPrintBase {
public:
    StringMetadataPrettyPrint( const Metadata& m ) : metadata_( m ) {}
    void print( std::ostream& out ) const override {
        std::string type = metadata_.getString( "type" );
        ATLAS_ASSERT( type == "string" );
        std::string value = metadata_.getString( "value" );
        if ( value.size() <= 32 ) {
            out << value;
        }
        else {
            out << value.substr( 0, 32 ) << "...";
        }
    }

private:
    Metadata metadata_;
};

class ScalarMetadataPrettyPrint : public MetadataPrettyPrintBase {
public:
    ScalarMetadataPrettyPrint( const Metadata& m ) : metadata_( m ) {}
    template <typename T>
    T decode() const {
        T value;
        Data dummy;
        atlas::io::decode( metadata_, dummy, value );
        return value;
    }
    void print( std::ostream& out ) const override {
        std::string type = metadata_.getString( "type" );
        ATLAS_ASSERT( type == "scalar" );
        std::string datatype = metadata_.getString( "datatype" );
        std::string base64   = metadata_.getString( "base64" );
        out << std::setw( 7 ) << std::left << datatype << ": ";
        if ( datatype == array::DataType::str<double>() ) {
            out << decode<double>();
        }
        else if ( datatype == array::DataType::str<float>() ) {
            out << decode<float>();
        }
        else if ( datatype == array::DataType::str<size_t>() ) {
            out << decode<size_t>();
        }
        else if ( datatype == array::DataType::str<std::int32_t>() ) {
            out << decode<std::int32_t>();
        }
        else if ( datatype == array::DataType::str<std::int64_t>() ) {
            out << decode<std::int64_t>();
        }
    }

private:
    Metadata metadata_;
};


MetadataPrettyPrint::MetadataPrettyPrint( const atlas::io::Metadata& m ) {
    std::string type = m.getString( "type" );

    // Poor man factory builder
    if ( type == "array" ) {
        impl_.reset( new ArrayMetadataPrettyPrint( m ) );
    }
    else if ( type == "scalar" ) {
        impl_.reset( new ScalarMetadataPrettyPrint( m ) );
    }
    else if ( type == "string" ) {
        impl_.reset( new StringMetadataPrettyPrint( m ) );
    }
    else {
        impl_.reset( new DefaultMetadataPrettyPrint( m ) );
    }
}

std::ostream& operator<<( std::ostream& out, const MetadataPrettyPrint& p ) {
    out << *p.impl_;
    return out;
}
std::string MetadataPrettyPrint::str() const {
    return impl_->str();
}


struct TablePrinter {
    std::vector<std::vector<std::string> > columns;
    std::vector<size_t> widths;
    std::string indent{"    "};
    std::string sep{" "};
    int col{0};
    int row{1};
    std::vector<bool> optional;
    std::vector<bool> underline;


    TablePrinter() = default;

    TablePrinter& column( const std::string& title, size_t width = 0 ) {
        columns.emplace_back( std::vector<std::string>{title} );
        widths.emplace_back( std::max( title.size(), width ) );
        optional.push_back( false );
        underline.push_back( true );
        return *this;
    }
    TablePrinter& column() {
        columns.emplace_back( std::vector<std::string>{""} );
        widths.emplace_back( 0 );
        optional.push_back( true );
        underline.push_back( false );
        return *this;
    }

    TablePrinter& operator<<( const MetadataPrettyPrint p ) {
        *this << p.str();
        return *this;
    }

    TablePrinter& operator<<( const std::string& s ) {
        columns[col].emplace_back( s );
        widths[col] = std::max( widths[col], s.size() );
        if ( optional[col] && widths[col] > 0 ) {
            optional[col] = false;
            widths[col]   = std::max( widths[col], columns[col][0].size() );
        }
        ++col;
        if ( col == columns.size() ) {
            ++row;
            col = 0;
        }
        return *this;
    }

    void print( std::ostream& out ) const {
        out << " ";
        for ( size_t c = 0; c < columns.size(); ++c ) {
            out << " " << sep << " " << std::setw( widths[c] ) << std::left << columns[c][0];
        }
        out << " " << sep << std::endl;

        out << " ";
        for ( size_t c = 0; c < columns.size(); ++c ) {
            const char u = underline[c] ? '-' : ' ';
            out << " " << sep << " " << std::string( widths[c], u );  // underline
        }
        out << " " << sep << std::endl;

        for ( size_t r = 1; r < row; ++r ) {
            out << " ";
            for ( size_t c = 0; c < columns.size(); ++c ) {
                out << " " << sep << " " << std::setw( widths[c] ) << std::left << columns[c][r];
            }
            out << " " << sep << std::endl;
        }
    }
};


TableFormat::TableFormat( const Record::URI& record, const util::Config& config ) :
    record_( Session::record( record.path, record.offset ) ) {
    for ( const auto& key : record_.keys() ) {
        items_.emplace( key, Metadata() );
        RecordItemReader{RecordItem::URI{record.path, record.offset, key}}.read( items_.at( key ) );
    }

    config.get( "details", print_details_ );
}

void TableFormat::print( std::ostream& out ) const {
    ATLAS_ASSERT( not record_.empty() );

    TablePrinter table;
    table.column( "name" );
    table.column( "type" );
    table.column( "description" );
    if ( print_details_ ) {
        table.column( "version" );
        table.column( "comp." );
        table.column( "size" );
        table.column( "checksum[:8]" );
        table.column( "endian" );
        table.column( "created" );
    }
    table.column();

    for ( auto& key : record_.keys() ) {
        auto& item          = items_.at( key );
        size_t cbytes       = item.data.compressed_size();
        size_t ubytes       = item.data.size();
        std::string endian  = ( item.data.endian() == Endian::little ) ? "little" : "big";
        std::string created = item.record.created().str().substr( 0, 16 );
        created[10]         = ' ';
        table << key;
        table << item.getString( "type" );
        table << MetadataPrettyPrint{item};
        if ( print_details_ ) {
            table << item.record.version();
            if ( item.data.section() ) {
                table << item.data.compression();
                table << Bytes{cbytes}.str() + ( item.data.compressed() ? " < " + Bytes{ubytes}.str() : "" );
                table << Checksum{item.data.checksum()}.str( 8 );
                table << endian;
            }
            else {
                table << ""
                      << ""
                      << ""
                      << "";
            }
            table << created;
        }

        if ( record_.metadata( key ).link() ) {
            table << "-> " + record_.metadata( key ).link().str();
        }
        else {
            table << "";
        }
    }
    table.print( out );
}


}  // namespace io
}  // namespace atlas
