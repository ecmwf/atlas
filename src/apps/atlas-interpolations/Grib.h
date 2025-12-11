#pragma once

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <memory>

#include "eccodes.h"

#include "eckit/filesystem/PathName.h"
#include "eckit/config/Configuration.h"


#include "atlas/grid.h"

// --------------------------------------------------------------------------------------------------------------

class Grib {
private:
    struct codes_handle_deleter {
        void operator()(codes_handle* h) {
            if (h) {
                codes_handle_delete(h);
            }
        }
    };

public:
    Grib(): grib_{nullptr, codes_handle_deleter()} {}

    Grib(codes_handle* h): grib_{h, codes_handle_deleter()} {}

    static Grib clone(const Grib& other) { return Grib(codes_handle_clone(other.handle())); }

    static Grib from_file(FILE* file) {
        int err;
        auto h = codes_handle_new_from_file(nullptr, file, PRODUCT_GRIB, &err);
        if (!h) {
            ATLAS_THROW_EXCEPTION("Could not create new grib handle from file. Error: " << err);
        }
        return Grib(h);
    }

    size_t get_values_size() const {
        size_t size;
        if (codes_get_size(handle(), "values", &size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get size of \"values\" from GRIB");
        };
        return size;
    }

    void get_values(double values[], size_t size) const {
        if (size < get_values_size()) {
            ATLAS_THROW_EXCEPTION("values[] has been allocated too small. encoded size: " << get_values_size() << ". allocated size: " << size);
        }
        if (codes_get_double_array(handle(), "values", values, &size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get \"values\" from GRIB");
        }
    }

    atlas::Grid get_grid() const {
        std::string gridType = get_string("gridType");
        if (gridType == "reduced_gg") {
            if (get_long("isOctahedral")) {
                return atlas::Grid("O" + std::to_string(get_long("N")));
            }
            else {
                return atlas::Grid("N" + std::to_string(get_long("N")));
            }
        }
        else if (gridType == "regular_gg") {
            return atlas::Grid("F" + std::to_string(get_long("N")));
        }
        else if (gridType == "regular_ll") {
            auto Ni   = get_long("Ni");
            auto Nj   = get_long("Nj");
            auto latN = get_double("latitudeOfFirstGridPointInDegrees");
            auto latS = get_double("latitudeOfLastGridPointInDegrees");
            auto lonW = get_double("longitudeOfFirstGridPointInDegrees");
            auto lonE = get_double("longitudeOfLastGridPointInDegrees");
            auto equal = [](double a, double b) -> bool {
                return std::abs(b-a) < 2.e-6;
            };

            if (equal(latN, 90.) && equal(latS, -90.) && equal(lonW, 0.) && equal(lonE, 360.-360./(double(Ni)))) {
                return atlas::Grid("L" + std::to_string(Ni) + "x" + std::to_string(Nj));
            }
            double lonShift = 360. / double(Ni) / 2.;
            double latShift = 180. / double(Nj) / 2.;
            if (equal(latN, 90.-latShift) && equal(latS, -90+latShift) && equal(lonW, lonShift) && equal(lonE, 360.-lonShift) ) {
                return atlas::Grid("S" + std::to_string(Ni) + "x" + std::to_string(Nj));
            }
            atlas::Log::error() << "Encoded regular_ll grid could not be validated as a global supported grid." << std::endl;
            atlas::Log::error() << "    - Ni:                                 " << Ni << std::endl;
            atlas::Log::error() << "    - Nj:                                 " << Nj << std::endl;
            atlas::Log::error() << "    - latitudeOfFirstGridPointInDegrees:  " << latN << std::endl;
            atlas::Log::error() << "    - latitudeOfLastGridPointInDegrees:   " << latS << std::endl;
            atlas::Log::error() << "    - longitudeOfFirstGridPointInDegrees: " << lonW << std::endl;
            atlas::Log::error() << "    - longitudeOfLastGridPointInDegrees:  " << lonE << std::endl;
            atlas::Log::error() << "    - iDirectionIncrementInDegrees:       " << get_double("iDirectionIncrementInDegrees") << std::endl;
            atlas::Log::error() << "    - jDirectionIncrementInDegrees:       " << get_double("jDirectionIncrementInDegrees") << std::endl;
        }
        else if(gridType == "healpix") {
            auto N = get_long("nside");
            auto orderingConvention = get_string("orderingConvention");
            if( orderingConvention == "ring" ) {
                return atlas::HealpixGrid(N);
            }
            atlas::Log::error() << "HealPIX ordering convention '" << orderingConvention << "' is not supported!" << std::endl;
        }
        return atlas::Grid();
    }

    std::string get_string(const std::string& key) const {
        char buffer[64];
        size_t size = sizeof(buffer);
        if (codes_get_string(handle(), key.c_str(), buffer, &size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get \"" << key << "\" from GRIB");
        }

        return std::string(buffer);
    };

    long get_long(const std::string& key) const {
        long value;
        if (codes_get_long(handle(), key.c_str(), &value) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get \"" << key << "\" from GRIB");
        }
        return value;
    };

    double get_double(const std::string& key) const {
        double value;
        if (codes_get_double(handle(), key.c_str(), &value) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get \"" << key << "\" from GRIB");
        }
        return value;
    };

    void set_grid(const atlas::Grid& agrid) {
        if( not atlas::StructuredGrid(agrid) ) {
            ATLAS_THROW_EXCEPTION("This grid is unsupported to write to GRIB: " << agrid.spec().json() );
        }

        auto grid = atlas::StructuredGrid(agrid);
        auto gaussian_grid = atlas::GaussianGrid(grid);
        auto regular_grid  = atlas::RegularGrid(grid);
        auto healpix_grid  = atlas::HealpixGrid(grid);

        set_string("gridType", atlas::RegularGaussianGrid(grid) ? "regular_gg"
                             : atlas::ReducedGaussianGrid(grid) ? "reduced_gg"
                             : atlas::RegularLonLatGrid(grid)   ? "regular_ll"
                             : atlas::HealpixGrid(grid)         ? "healpix"
                                                                : "unknown");

        if( healpix_grid ) {
            set_long("nside", healpix_grid.N() );
            set_string("orderingConvention","ring");
        }
        else {
            if (gaussian_grid) {
                set_long("N", gaussian_grid.N());
            }

            if (grid.reduced()) {
                set_long_array("pl", grid.nx());
            }

            if (regular_grid) {
                set_long("Ni", regular_grid.nx());
                set_long("Nj", regular_grid.ny());
                set_double("iDirectionIncrementInDegrees",360./double(regular_grid.nx()));
                if( atlas::RegularLonLatGrid() ) {
                    set_double("jDirectionIncrementInDegrees",std::abs(regular_grid.y(1)-regular_grid.y(0)));
                }
            }
            set_double("latitudeOfFirstGridPointInDegrees", grid.y(0));
            set_double("latitudeOfLastGridPointInDegrees", grid.y(grid.ny() - 1));
            set_double("longitudeOfFirstGridPointInDegrees", grid.x(0, 0));
            set_double("longitudeOfLastGridPointInDegrees", grid.x(0, 0) + 360. - 360. / double(grid.nxmax()));
        }
    }

    void set_values(const double values[], const size_t size) {
        if (codes_set_double_array(handle(), "values", values, size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not set \"values\" in GRIB");
        }
    }

    codes_handle* handle() const { return grib_.get(); }

    void set_string(const std::string& key, const std::string& value) {
        size_t size = value.size();
        if (codes_set_string(handle(), key.c_str(), value.c_str(), &size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not set {\"" << key << "\" : \"" << value << "\"} in GRIB");
        }
    }

    void set_long(const std::string& key, long value) {
        if (codes_set_long(handle(), key.c_str(), value) != 0) {
            ATLAS_THROW_EXCEPTION("Could not set \"" << key << "\" in GRIB");
        }
    }

    void set_double(const std::string& key, double value) {
        if (codes_set_double(handle(), key.c_str(), value) != 0) {
            ATLAS_THROW_EXCEPTION("Could not set \"" << key << "\" in GRIB");
        }
    }

    template <typename T>
    void set_long_array(const std::string& key, const std::vector<T>& value) {
        std::vector<long> v(value.begin(), value.end());
        if (codes_set_long_array(handle(), key.c_str(), v.data(), v.size()) != 0) {
            ATLAS_THROW_EXCEPTION("Could not set \"" << key << "\" in GRIB");
        }
    }


private:
    std::shared_ptr<codes_handle> grib_ = nullptr;
};

// --------------------------------------------------------------------------------------------------------------
