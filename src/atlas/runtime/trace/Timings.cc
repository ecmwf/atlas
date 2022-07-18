/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "Timings.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <regex>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/trace/CallStack.h"
#include "atlas/runtime/trace/CodeLocation.h"
#include "atlas/util/Config.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class TimingsRegistry {
private:
    std::vector<long> counts_;
    std::vector<double> tot_timings_;
    std::vector<double> min_timings_;
    std::vector<double> max_timings_;
    std::vector<double> var_timings_;
    std::vector<std::string> titles_;
    std::vector<CodeLocation> locations_;
    std::vector<long> nest_;
    std::vector<CallStack> stack_;
    std::map<size_t, size_t> index_;

    std::map<std::string, std::vector<size_t>> labels_;

    TimingsRegistry() = default;

public:
    static TimingsRegistry& instance() {
        static TimingsRegistry registry;
        return registry;
    }

    size_t add(const CodeLocation&, const CallStack& stack, const std::string& title, const Timings::Labels&);

    void update(size_t idx, double seconds);

    size_t size() const;

    void report(std::ostream& out, const eckit::Configuration& config);

private:
    std::string filter_filepath(const std::string& filepath) const;

    friend class Tree;
    friend class Node;
};

struct Node {
    Node(): index(-1) {}
    Node(size_t _index): index(_index) {
        size_t _nest = TimingsRegistry::instance().nest_[index];


        const auto& this_stack = TimingsRegistry::instance().stack_[index];
        auto this_stack_hash   = TimingsRegistry::instance().stack_[index].hash();
        auto is_child          = [&](size_t i) -> bool {
            CallStack child_stack = TimingsRegistry::instance().stack_[i];
            child_stack.pop();
            auto child_stack_hash = child_stack.hash();
            return child_stack_hash == this_stack_hash;
        };


        for (size_t i = index + 1; i < TimingsRegistry::instance().size(); ++i) {
            if (TimingsRegistry::instance().nest_[i] == _nest + 1) {
                if (is_child(i)) {
                    children.emplace_back(new Node(i));
                }
            }
        }
    }
    std::vector<std::unique_ptr<Node>> children;
    std::unique_ptr<Node> parent;
    long index;
    void print(std::ostream& out) {
        if (index >= 0) {
            size_t _nest = nest();
            for (size_t i = 1; i < _nest; ++i) {
                out << "    ";
            }
            out << TimingsRegistry::instance().titles_[index] << std::endl;
        }
        for (auto& child : children) {
            child->print(out);
        }
    }
    size_t nest() const { return TimingsRegistry::instance().nest_[index]; }
    void order(std::vector<size_t>& order) const {
        if (index >= 0) {
            order.emplace_back(index);
        }
        for (auto& child : children) {
            child->order(order);
        }
    }
};
struct Tree {
    Node root;
    Tree() {
        for (size_t j = 0; j < TimingsRegistry::instance().size(); ++j) {
            auto& nest = TimingsRegistry::instance().nest_[j];
            if (nest == 1) {
                auto& children = root.children;
                children.emplace_back(new Node(j));
            }
        }
    }
    void print(std::ostream& out) { root.print(out); }
    std::vector<size_t> order() const {
        std::vector<size_t> order;
        order.reserve(TimingsRegistry::instance().size());
        root.order(order);
        ATLAS_ASSERT(order.size() == TimingsRegistry::instance().size(),
                     "Likely a atlas_Trace has not finalised properly");
        return order;
    }
};

size_t TimingsRegistry::add(const CodeLocation& loc, const CallStack& stack, const std::string& title,
                            const Timings::Labels& labels) {
    size_t key = stack.hash();
    auto it    = index_.find(key);
    if (it == index_.end()) {
        size_t idx  = size();
        index_[key] = idx;
        counts_.emplace_back(0);
        tot_timings_.emplace_back(0);
        min_timings_.emplace_back(std::numeric_limits<double>::max());
        max_timings_.emplace_back(0);
        var_timings_.emplace_back(0);
        titles_.emplace_back(title);
        locations_.emplace_back(loc);
        nest_.emplace_back(stack.size());
        stack_.emplace_back(stack);

        for (const auto& label : labels) {
            labels_[label].emplace_back(idx);
        }

        return idx;
    }
    else {
        return it->second;
    }
}

void TimingsRegistry::update(size_t idx, double seconds) {
    auto sqr          = [](double x) { return x * x; };
    double n          = counts_[idx] + 1;
    double avg_nm1    = tot_timings_[idx] / std::max(n, 1.);
    double var_nm1    = var_timings_[idx];
    var_timings_[idx] = n == 1. ? 0. : (n - 2.) / (n - 1.) * var_nm1 + 1. / n * sqr(seconds - avg_nm1);
    min_timings_[idx] = std::min(seconds, min_timings_[idx]);
    max_timings_[idx] = std::max(seconds, max_timings_[idx]);
    tot_timings_[idx] += seconds;
    counts_[idx] += 1;
}

size_t TimingsRegistry::size() const {
    return counts_.size();
}

void TimingsRegistry::report(std::ostream& out, const eckit::Configuration& config) {
    auto box_horizontal = [](int n) {
        std::string s;
        s.reserve(2 * n);
        for (int i = 0; i < n; ++i) {
            s += "\u2500";
        }
        return s;
    };
    std::string box_corner_tl("\u250c");
    std::string box_corner_tr("\u2510");
    std::string box_corner_bl("\u2514");
    std::string box_corner_br("\u2518");
    std::string box_vertical("\u2502");
    std::string box_T_down("\u252C");
    std::string box_T_up("\u2534");
    std::string box_T_right("\u251C");
    std::string box_T_left("\u2524");
    std::string box_cross("\u253C");

    long indent                                     = config.getLong("indent", 2);
    long depth                                      = config.getLong("depth", 0);
    long decimals                                   = config.getLong("decimals", 5);
    bool header                                     = config.getBool("header", true);
    std::vector<std::string> excluded_labels_vector = config.getStringVector("exclude", std::vector<std::string>());
    std::vector<std::string> include_back;

    auto order = Tree().order();

    for (auto& label : excluded_labels_vector) {
        size_t found = label.find("/*");
        if (found != std::string::npos) {
            label.erase(found, 2);
            include_back.push_back(label);
        }
    }

    std::set<std::string> excluded_labels(excluded_labels_vector.begin(), excluded_labels_vector.end());

    auto digits_before_decimal = [](double x) -> int {
        return std::floor(std::log10(std::trunc(std::max(1., x)))) + 1;
    };
    auto digits = [](long x) -> long { return std::floor(std::log10(std::max(1l, x))) + 1l; };

    std::vector<size_t> excluded_timers_vector;
    for (auto& label : labels_) {
        auto name = label.first;
        if (excluded_labels.count(name)) {
            auto timers = label.second;
            for (size_t j : timers) {
                excluded_timers_vector.push_back(j);
            }
        }
    }
    std::set<size_t> excluded_timers(excluded_timers_vector.begin(), excluded_timers_vector.end());

    auto excluded = [&](size_t i) -> bool {
        if (depth and nest_[i] > depth) {
            return true;
        }
        return excluded_timers.count(i);
    };

    std::vector<long> excluded_nest_stored(size());
    long excluded_nest = size();
    for (size_t jj = 0; jj < size(); ++jj) {
        size_t j = order[jj];
        if (nest_[j] > excluded_nest) {
            excluded_timers.insert(j);
        }
        if (not excluded(j)) {
            excluded_nest = nest_[j] + 1;
        }
        else {
            excluded_nest = std::min(excluded_nest, nest_[j]);
        }
        excluded_nest_stored[j] = excluded_nest;
    }
    for (auto& label : include_back) {
        auto timers = labels_[label];
        for (size_t j : timers) {
            if (nest_[j] == excluded_nest_stored[j]) {
                excluded_timers.erase(j);
            }
        }
    }

    size_t max_title_length(0);
    size_t max_location_length(0);
    size_t max_nest(0);
    long max_count(0);
    double max_seconds(0);
    for (size_t j = 0; j < size(); ++j) {
        size_t nest = nest_[j];
        max_nest    = std::max(max_nest, nest);
        if (not excluded(j)) {
            const auto& loc        = locations_[j];
            max_title_length       = std::max(max_title_length, titles_[j].size() + nest_[j] * indent);
            max_count              = std::max(max_count, counts_[j]);
            max_seconds            = std::max(max_seconds, tot_timings_[j]);
            size_t location_length = filter_filepath(loc.file()).size() + 2 + digits(loc.line());
            max_location_length    = std::max(max_location_length, location_length);
        }
    }
    size_t max_count_length = digits(max_count);
    if (header) {
        max_count_length = std::max(std::string("cnt").size(), max_count_length);
    }
    size_t max_digits_before_decimal = digits_before_decimal(max_seconds);

    auto print_time = [max_digits_before_decimal, decimals](double x) -> std::string {
        std::stringstream out;
        char unit = 's';
        if (std::floor(x) >= 60) {
            x /= 60.;
            unit = 'm';
        }
        out << std::right << std::fixed << std::setprecision(decimals)
            << std::setw(max_digits_before_decimal + decimals + 1) << x << unit;
        return out.str();
    };

    auto print_line = [&](size_t length) -> std::string { return box_horizontal(length); };

    auto print_horizontal = [&](const std::string& sep) -> std::string {
        std::stringstream ss;
        ss << print_line(max_title_length + digits(size()) + 3) << sep << print_line(max_count_length) << sep
           << print_line(max_digits_before_decimal + decimals + 2) << sep
           << print_line(max_digits_before_decimal + decimals + 2) << sep
           << print_line(max_digits_before_decimal + decimals + 2) << sep
           << print_line(max_digits_before_decimal + decimals + 2) << sep
           << print_line(max_digits_before_decimal + decimals + 2) << sep << print_line(max_location_length);
        return ss.str();
    };

    std::string sept = box_horizontal(1) + box_T_down + box_horizontal(1);
    std::string seph = box_horizontal(1) + box_cross + box_horizontal(1);
    std::string sep  = std::string(" ") + box_vertical + std::string(" ");
    std::string sepf = box_horizontal(1) + box_T_up + box_horizontal(1);

    out << print_horizontal(sept) << std::endl;
    out << std::left << std::setw(max_title_length + digits(size()) + 3) << "Timers" << sep
        << std::setw(max_count_length) << "cnt" << sep << std::setw(max_digits_before_decimal + decimals + 2ul) << "tot"
        << sep << std::setw(max_digits_before_decimal + decimals + 2ul) << "avg" << sep
        << std::setw(max_digits_before_decimal + decimals + 2ul) << "std" << sep
        << std::setw(max_digits_before_decimal + decimals + 2ul) << "min" << sep
        << std::setw(max_digits_before_decimal + decimals + 2ul) << "max" << sep << "location" << std::endl;
    out << print_horizontal(seph) << std::endl;

    std::vector<std::string> prefix_(size());
    if (indent) {
        std::vector<bool> active(max_nest, false);
        for (long kk = long(size()) - 1; kk >= 0; --kk) {
            long k           = order[kk];
            const auto& nest = nest_[k];

            const CallStack& this_stack = stack_[k];
            const CallStack* next_stack_ptr;
            if (kk == size() - 1) {
                next_stack_ptr = &this_stack;
            }
            else {
                next_stack_ptr = &stack_[order[kk + 1]];
            }
            const CallStack& next_stack = *next_stack_ptr;

            auto this_it = this_stack.begin();
            auto next_it = next_stack.begin();
            for (size_t i = 0; this_it != this_stack.end() && next_it != next_stack.end(); ++i, ++this_it, ++next_it) {
                if (*this_it == *next_it) {
                    active[i] = active[i] or false;
                }
                else {
                    active[i] = true;
                }
            }
            for (size_t i = nest; i < active.size(); ++i) {
                active[i] = false;
            }

            std::stringstream out;
            for (long i = 0; i < nest - 1; ++i) {
                if (active[i]) {
                    out << box_vertical;
                }
                else {
                    out << " ";
                }
                for (long j = 1; j < indent; ++j) {
                    out << " ";
                }
            }
            if (active[nest - 1]) {
                out << box_T_right;
            }
            else {
                out << box_corner_bl;
            }
            for (long j = 1; j < indent; ++j) {
                out << box_horizontal(1);
            }

            prefix_[k] = out.str();
        }
    }

    for (size_t i = 0; i < size(); ++i) {
        size_t j    = order[i];
        auto& tot   = tot_timings_[j];
        auto& max   = max_timings_[j];
        auto& min   = std::min(max, min_timings_[j]);
        auto& count = counts_[j];
        auto& title = titles_[j];
        auto& loc   = locations_[j];
        auto& nest  = nest_[j];
        auto std    = std::sqrt(var_timings_[j]);
        auto avg    = (count == 0 ? 0. : tot / double(count));

        // mpi::comm().allReduceInPlace(min,eckit::mpi::min());
        // mpi::comm().allReduceInPlace(max,eckit::mpi::max());

        if (not excluded(j)) {
            out << std::setw(digits(long(size()))) << j << " : " << prefix_[j]  // prefix(indent,nest,next_nest)
                << std::left << std::setw(max_title_length - nest * indent) << title << sep
                << std::string(header ? "" : "count: ") << std::left << std::setw(max_count_length) << count << sep
                << std::string(header ? "" : "tot: ") << print_time(tot) << sep << std::string(header ? "" : "avg: ")
                << print_time(avg) << sep << std::string(header ? "" : "std: ") << print_time(std) << sep
                << std::string(header ? "" : "min: ") << print_time(min) << sep << std::string(header ? "" : "max: ")
                << print_time(max) << sep << filter_filepath(loc.file()) << " +" << loc.line() << std::endl;
        }
    }

    out << print_horizontal(sepf) << std::endl;

    std::string sepc = box_horizontal(3);

    out << std::left << box_horizontal(40) << sept << box_horizontal(5) << sept << box_horizontal(12) << "\n";
    out << std::left << std::setw(40) << "Timers accumulated by label" << sep << std::left << std::setw(5) << "count"
        << sep << "time" << std::endl;
    out << std::left << box_horizontal(40) << seph << box_horizontal(5) << seph << box_horizontal(12) << "\n";
    for (auto& label : labels_) {
        auto name   = label.first;
        auto timers = label.second;
        double tot(0);
        double count(0);
        for (size_t j : timers) {
            tot += tot_timings_[j];
            count += counts_[j];
        }
        out << std::left << std::setw(40) << name << sep << std::left << std::setw(5) << count << sep << print_time(tot)
            << std::endl;
    }
    out << std::left << box_horizontal(40) << sepf << box_horizontal(5) << sepf << box_horizontal(12) << "\n";
}

std::string TimingsRegistry::filter_filepath(const std::string& filepath) const {
    std::smatch matches;
    std::string basename = eckit::PathName(filepath).baseName();
    if (std::regex_search(filepath, matches, std::regex{"(.*)?/atlas/src/(.*)"})) {
        return "[atlas] " + basename;
    }
    if (std::regex_search(filepath, matches, std::regex{"(.*)?/atlas-io/src/(.*)"})) {
        return "[atlas-io] " + basename;
    }
    return basename;
}

Timings::Identifier Timings::add(const CodeLocation& loc, const CallStack& stack, const std::string& title,
                                 const Labels& labels) {
    return TimingsRegistry::instance().add(loc, stack, title, labels);
}

void Timings::update(const Identifier& id, double seconds) {
    TimingsRegistry::instance().update(id, seconds);
}

std::string Timings::report() {
    return report(util::NoConfig());
}

std::string Timings::report(const Configuration& config) {
    std::ostringstream out;
    TimingsRegistry::instance().report(out, config);
    return out.str();
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
