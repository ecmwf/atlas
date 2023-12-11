/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "GridPointsJSONWriter.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "eckit/utils/Tokenizer.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

namespace {
std::vector<long> points_from_list(const std::string& list, long base) {
    std::vector<long> points;
    if (list.empty()) {
        return points;
    }
    if (list[0] == '[' && list[list.size()-1] == ']') {
        return points_from_list(std::string(list.begin()+1,list.begin()+list.size()-1), base);
    }
    std::vector<std::string> points_ranges;
    eckit::Tokenizer{","}(list,points_ranges);
    auto tokenize_range = eckit::Tokenizer{"-"};
    auto to_int = [](const std::string& s, const eckit::CodeLocation& here) -> long {
        try {
            return std::stol(s);
        }
        catch( std::exception& e) {
            throw_Exception("Could not convert '" + s + "' to integer",here);
        }
    };
    for (const auto& point_range: points_ranges) {
        std::vector<std::string> point_range_split;
        tokenize_range(point_range,point_range_split);
        if (point_range_split.size() == 1) {
            points.emplace_back(to_int(point_range,Here())-base);
        }
        else {
            for (long p = to_int(point_range_split[0],Here()); p <= to_int(point_range_split[1],Here()); ++p) {
                points.emplace_back(p-base);
            }
        }
    }
    return points;
}
}

//------------------------------------------------------------------------------

GridPointsJSONWriter::GridPointsJSONWriter(Grid grid, const eckit::Parametrisation& args) : grid_{grid} {
    args.get("json.precision",precision_=-1);
    args.get("verbose",verbose_=0);
    if (not args.get("partitions",nb_partitions_=0)) {
        args.get("partitioner.partitions",nb_partitions_);
    }
    if (not args.get("partition",partition_=-1)) {
        args.get("partition",partition_);
    }
    args.get("json.pretty", pretty_=false);
    args.get("field",field_="lonlat");
    args.get("field_base",field_base_=0);
    std::string points_list;
    if (args.get("index",points_list)) {
        args.get("index_base",points_base_ = 0);
        points_ = points_from_list(points_list,points_base_);
    }

    if( nb_partitions_ > 0 ) {
        auto partitioner_config = grid_.partitioner();  // recommended by the grid itself
        std::string partitioner;
        if (args.get("partitioner.type",partitioner)) {
            partitioner_config.set("type", partitioner);
        }
        partitioner_config.set("partitions", nb_partitions_);
        distribution_ = grid::Distribution{grid_,partitioner_config};
    }
}

//------------------------------------------------------------------------------

void GridPointsJSONWriter::write(std::ostream& out, eckit::Channel& info) const {
    write(out, &info);
}

//------------------------------------------------------------------------------

void GridPointsJSONWriter::write(std::ostream& out, std::ostream* info) const {

    if (field_ == "none") {
        return;
    }

    int points_newline = pretty_ ? true : false;
    int points_indent = pretty_ ? 2 : 0;
    int partition_indent = 0;
    size_t chunk_size = 1000000;

    auto saved_fmtflags = out.flags();
    if (precision_ >= 0) {
        out << std::fixed << std::setprecision(precision_);
    }

    auto writePoint = [&out](const Point2& p) {
        out << "[" << p[0] << "," << p[1] << "]";
    };

    if(nb_partitions_ == 0 || field_ == "partition") {
        out << "[";
        if (pretty_) {
            out << '\n';
        }
        std::string join = points_newline ? ",\n" : ",";
        std::string indent(points_indent,' ');

        size_t n{0};
        if (points_.size()) {
            if (field_ == "lonlat") {
                auto lonlat = grid_.lonlat();
                for (auto p: points_) {
                    if (p < grid_.size()) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        auto it = lonlat.begin()+(p);
                        writePoint(*it);
                        ++n;
                    }
                }
            }
            else if (field_ == "index") {
                for (auto p: points_) {
                    if (p < grid_.size()) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        out << p + field_base_;
                        ++n;
                    }
                }
            }
            else if (field_ == "partition") {
                for (auto p: points_) {
                    if (p < grid_.size()) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        out << distribution_.partition(p) + field_base_;
                        ++n;
                    }
                }
            }
            else {
                ATLAS_THROW_EXCEPTION("Cannot output field \""<<field_<<"\"");
            }
        }
        else {
            size_t size = grid_.size();
            size_t end = 0;
            size_t begin;
            auto write_progress = [&]() {
                if (info && verbose_ && grid_.size() > chunk_size) {
                    *info << std::fixed << std::setprecision(1) << double(begin)/double(grid_.size()) * 100 << "% completed" << std::endl;
                }
            };
            auto write_progress_end = [&]() {
                if (info && verbose_ && grid_.size() > chunk_size) {
                    *info << "100% completed" << std::endl;
                }
            };
            while(end != grid_.size()) {
                begin = end;
                end = std::min<size_t>(grid_.size(),begin+chunk_size);
                write_progress();
                if (field_ == "lonlat") {
                    auto it = grid_.lonlat().begin()+begin;
                    for (size_t n=begin; n<end; ++n, ++it) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        writePoint(*it);
                    }
                }
                else if (field_ == "index" ) {
                    for (size_t n=begin; n<end; ++n) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        out << n + field_base_;
                    }
                }
                else if (field_ == "partition" ) {
                    for (size_t n=begin; n<end; ++n) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        out << distribution_.partition(n) + field_base_;
                    }
                }
                else {
                    ATLAS_THROW_EXCEPTION("Cannot output field \""<<field_<<"\"");
                }
            }
            write_progress_end();
        }

        if (pretty_) {
            out << '\n';
        }
        out << "]";
        out << std::endl;
    }
    else {
        std::vector<int> part(std::min<size_t>(chunk_size,grid_.size()));
        auto write_chunk = [&](int partition, size_t begin, size_t end, size_t& n) {
            size_t size{end - begin};
            if (partition != -1) {
                distribution_.partition(begin, end, part);
                auto part_begin = part.begin();
                auto part_end = part.begin()+size;
                if (std::find(part_begin,part_end,partition) == part_end) {
                    return; // no points in this chunk to write
                }
            }
            else {
                part.assign(part.size(),-1);
            }
            std::string join = points_newline ? ",\n" : ",";
            std::string indent(points_indent,' ');
            if( field_ == "lonlat" ) {
                auto it = grid_.lonlat().begin() + begin;
                for (size_t j=0; j<size; ++j, ++it) {
                    if (part[j] == partition) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        writePoint(*it);
                        ++n;
                    }
                }
            }
            else if (field_ == "index") {
                for (size_t j=0; j<size; ++j) {
                    if (part[j] == partition) {
                        if (n!=0) {
                            out << join;
                        }
                        out << indent;
                        out << begin + j + field_base_;
                        ++n;
                    }
                }
            }
            else {
                ATLAS_THROW_EXCEPTION("Cannot output field \""<<field_<<"\"");
            }
        };

        auto write_partition = [&](int partition) {
            if (info && verbose_) {
                *info << "Partition " << partition << std::endl;
            }
            size_t end = 0;
            size_t n{0};
            out << std::string(partition_indent,' ') << "[";
            if (pretty_) {
                out << '\n';
            }
            while(end != grid_.size()) {
                auto begin = end;
                if (info && verbose_ && grid_.size() > chunk_size) {
                    *info << std::fixed << std::setprecision(1) << double(begin)/double(grid_.size()) * 100 << "% completed" << std::endl;
                }
                end = std::min<size_t>(grid_.size(),begin+chunk_size);
                write_chunk(partition,begin,end,n);
            }
            if (verbose_) {
                if (info && grid_.size() > chunk_size) {
                    *info << "100% completed" << std::endl;
                }
            }
            if (pretty_) {
                out << '\n';
            }
            out << std::string(partition_indent,' ') << "]" << std::flush;
            return n;
        };

        if (partition_ >= 0) {
            auto n = write_partition(partition_);
            out << std::endl;
            if (info) {
                *info << "Partition " << partition_ << " contains " << n << " points." << std::endl;
            }
        }
        else {
            if (pretty_) {
                partition_indent += 2;
                points_indent += 2;
            }
            out << "[\n";
            std::vector<size_t> points_per_partition(nb_partitions_);
            for( int p = 0; p < nb_partitions_; ++p ) {
                if (p != 0) {
                    out << ",\n";
                }
                points_per_partition[p] = write_partition(p);
            }
            out << "\n]";
            out << std::endl;
            if (info) {
                *info << "Points per partition:" << std::endl;
                for (size_t p=0; p<nb_partitions_; ++p) {
                    *info << "    " << std::setw(5) << std::left << p << points_per_partition[p] << std::endl;
                }
            }
        }
    }
    out.flags(saved_fmtflags);
}

//------------------------------------------------------------------------------

}
}