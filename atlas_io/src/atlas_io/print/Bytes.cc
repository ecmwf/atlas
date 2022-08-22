/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Bytes.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace atlas {
namespace io {


namespace {

template <typename Value>
class FixedFormat {
public:
    using value_type = Value;
    FixedFormat(value_type x, long precision): x_(x), precision_(precision > 0 ? precision : 20) {}
    void print(std::ostream& out) const {
        for (long precision = 0; precision <= precision_; ++precision) {
            if (is_precision(precision) || precision == precision_) {
                out << std::setprecision(precision);
                out << std::fixed << x_;
                break;
            }
        }
    }

    bool is_precision(long precision) const {
        std::stringstream ss;
        ss << std::setprecision(precision);
        ss << std::fixed << x_;
        value_type _x;
        ss >> _x;
        return std::abs(x_ - _x) < 1.e-20;
    }

    friend std::ostream& operator<<(std::ostream& out, const FixedFormat& This) {
        This.print(out);
        return out;
    }

private:
    value_type x_;
    long precision_;
};

inline FixedFormat<double> fixed_format(double x, long precision) {
    return FixedFormat<double>(x, precision);
}

/*
static std::pair<double, std::string> reduce_to_10( size_t bytes ) {
    static const std::vector<std::string> magnitudes{"B", "K", "M", "G", "T", "P", "E", "Z", "Y"};
    double x = bytes;
    size_t n = 0;
    while ( x >= 10 && n < magnitudes.size() ) {
        x /= 1024.;
        n++;
    }
    return std::make_pair( x, magnitudes[n] );
}
*/

static std::pair<double, std::string> reduce_to_1000(size_t bytes) {
    static const std::vector<std::string> magnitudes{"B", "K", "M", "G", "T", "P", "E", "Z", "Y"};
    double x = bytes;
    size_t n = 0;
    while (x >= 1000 && n < magnitudes.size()) {
        x /= 1024.;
        n++;
    }
    return std::make_pair(x, magnitudes[n]);
}


}  // namespace

void Bytes::print(std::ostream& out, int decimals, int width) const {
    if (bytes_ < 1000 && width >= 4) {
        out << std::setw(width - 1) << std::right << bytes_ << 'B';
    }
    else {
        auto pair = reduce_to_1000(bytes_);
        out << std::setw(width - 1) << std::right << fixed_format(pair.first, decimals);
        out << pair.second;
    }
}

std::ostream& operator<<(std::ostream& out, const Bytes& bytes) {
    bytes.print(out);
    return out;
}


std::string Bytes::str(int decimals, int width) const {
    std::stringstream s;
    print(s, decimals, width);
    return s.str();
}


}  // namespace io
}  // namespace atlas
