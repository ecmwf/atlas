/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck

#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {

class GridBuilder {
public:
    using Registry = std::map<std::string, GridBuilder*>;
    static const Registry& nameRegistry();
    static const Registry& typeRegistry();

public:
    GridBuilder(const std::string& names);
    GridBuilder(const std::string& type, const std::vector<std::string>& regexes,
                const std::vector<std::string>& names);

    virtual ~GridBuilder();

    virtual const Grid::Implementation* create(const Grid::Config&) const;

    virtual const Grid::Implementation* create(const std::string&, const Grid::Config& = Grid::Config()) const = 0;

    virtual const std::string& type() const;
    const std::vector<std::string>& regexes() const;
    const std::vector<std::string>& names() const;

    void registerNamedGrid(const std::string& name);

protected:
    bool match(const std::string& string, std::vector<std::string>& matches, int& id) const;

private:
    friend std::ostream& operator<<(std::ostream& os, const GridBuilder& g);

    virtual void print(std::ostream& os) const = 0;

    std::vector<std::string> names_;
    std::vector<std::string> pretty_names_;
    std::string type_;
};

}  // namespace grid
}  // namespace atlas
