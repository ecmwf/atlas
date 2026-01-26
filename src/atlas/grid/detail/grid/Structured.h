/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <memory>

#include "atlas/grid/Spacing.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

/**
 * @brief Structured Grid
 *
 * This class is a base class for all grids that can be described by
 * constant latitudes with a uniform distribution of points per latitude
 * in zonal direction.
 * This means any full grid and reduced grid, both regular, gaussian or other
 * such distribution can be represented with this class
 */
class Structured : public Grid {
private:
    struct ComputePointXY {
        ComputePointXY(const Structured& grid): grid_(grid), ny_(grid_.ny()) {}
        void operator()(idx_t i, idx_t j, PointXY& point) {
            if (j < ny_) {  // likely
                grid_.xy(i, j, point.data());
            }
        }
        const Structured& grid_;
        idx_t ny_;
    };
    struct ComputePointLonLat {
        ComputePointLonLat(const Structured& grid): grid_(grid), ny_(grid_.ny()) {}
        void operator()(idx_t i, idx_t j, PointLonLat& point) {
            if (j < ny_) {  // likely
                grid_.lonlat(i, j, point.data());
            }
        }
        const Structured& grid_;
        idx_t ny_;
    };

    template <typename Base, typename ComputePoint>
    class StructuredIterator : public Base {
    public:
        StructuredIterator(const Structured& grid, bool begin = true):
            grid_(grid), ny_(grid_.ny()), i_(0), j_(begin ? 0 : grid_.ny()), compute_point{grid_} {
            if (j_ != ny_ && grid_.size()) {
                compute_point(i_, j_, point_);
            }
        }
        virtual bool next(typename Base::value_type& point) {
            if (j_ < ny_ && i_ < grid_.nx(j_)) {
                compute_point(i_++, j_, point);

                if (i_ == grid_.nx(j_)) {
                    j_++;
                    i_ = 0;
                }
                return true;
            }
            return false;
        }

        virtual const typename Base::reference operator*() const { return point_; }

        virtual const Base& operator++() {
            ++i_;
            if (i_ == grid_.nx(j_)) {
                ++j_;
                i_ = 0;
            }
            compute_point(i_, j_, point_);
            return *this;
        }

        virtual const Base& operator+=(typename Base::difference_type distance) {
            idx_t d = distance;
            while (j_ != ny_ && d >= (grid_.nx(j_) - i_)) {
                d -= (grid_.nx(j_) - i_);
                ++j_;
                i_ = 0;
            }
            i_ += d;
            compute_point(i_, j_, point_);
            return *this;
        }

        virtual typename Base::difference_type distance(const Base& other) const {
            const auto& _other               = static_cast<const StructuredIterator&>(other);
            typename Base::difference_type d = 0;
            idx_t j                          = j_;
            idx_t i                          = i_;
            while (j < _other.j_) {
                d += grid_.nx(j) - i;
                ++j;
                i = 0;
            }
            d += _other.i_;
            return d;
        }

        virtual bool operator==(const Base& other) const {
            return j_ == static_cast<const StructuredIterator&>(other).j_ &&
                   i_ == static_cast<const StructuredIterator&>(other).i_;
        }

        virtual bool operator!=(const Base& other) const {
            return i_ != static_cast<const StructuredIterator&>(other).i_ ||
                   j_ != static_cast<const StructuredIterator&>(other).j_;
        }

        virtual std::unique_ptr<Base> clone() const {
            auto result    = new StructuredIterator(grid_, false);
            result->i_     = i_;
            result->j_     = j_;
            result->point_ = point_;
            return std::unique_ptr<Base>(result);
        }

    public:
        const Structured& grid_;
        idx_t ny_;
        idx_t i_;
        idx_t j_;
        typename Base::value_type point_;
        ComputePoint compute_point;
    };

public:
    using IteratorXY     = StructuredIterator<Grid::IteratorXY, ComputePointXY>;
    using IteratorLonLat = StructuredIterator<Grid::IteratorLonLat, ComputePointLonLat>;

public:
    class XSpace {
        class Implementation : public util::Object {
        public:
            // Constructor NVector can be either std::vector<int> or std::vector<long>
            template <typename NVector>
            Implementation(const std::array<double, 2>& interval, const NVector& N, bool endpoint = true);

            Implementation(const std::array<double, 2>& interval, std::initializer_list<int>&& N, bool endpoint = true);

            Implementation(const Spacing&, idx_t ny = 1);

            Implementation(const std::vector<Spacing>&);

            Implementation(const Config&);

            Implementation(const std::vector<Config>&);

            Implementation(const std::vector<spacing::LinearSpacing::Params>&);

            idx_t ny() const { return ny_; }

            // Minimum number of points across parallels (constant y)
            idx_t nxmin() const { return nxmin_; }

            // Maximum number of points across parallels (constant y)
            idx_t nxmax() const { return nxmax_; }

            /// Number of points per latitude
            const std::vector<idx_t>& nx() const { return nx_; }

            /// Value of minimum longitude per latitude [default=0]
            const std::vector<double>& xmin() const { return xmin_; }

            /// Value of maximum longitude per latitude [default=0]
            const std::vector<double>& xmax() const { return xmax_; }

            /// Value of longitude increment
            const std::vector<double>& dx() const { return dx_; }

            /// Value of minimum x over entire grid
            double min() const { return min_; }

            /// Value of maximum x over entire grid
            double max() const { return max_; }

            Spec spec() const;

            std::string type() const;

        private:
            void reserve(idx_t ny);

        private:
            idx_t ny_;
            idx_t nxmin_;
            idx_t nxmax_;
            std::vector<idx_t> nx_;
            std::vector<double> xmin_;
            std::vector<double> xmax_;
            std::vector<double> dx_;
            double min_;
            double max_;
        };

    public:
        XSpace();

        XSpace(const Spacing&, idx_t ny = 1);

        XSpace(const std::vector<Spacing>&);

        XSpace(const std::vector<spacing::LinearSpacing::Params>&);

        // Constructor NVector can be either std::vector<int> or std::vector<long> or initializer list
        template <typename NVector>
        XSpace(const std::array<double, 2>& interval, const NVector& N, bool endpoint = true);

        XSpace(const std::array<double, 2>& interval, std::initializer_list<int>&& N, bool endpoint = true);

        XSpace(const Config&);

        XSpace(const std::vector<Config>&);

        idx_t ny() const { return impl_->ny(); }

        // Minimum number of points across parallels (constant y)
        idx_t nxmin() const { return impl_->nxmin(); }

        // Maximum number of points across parallels (constant y)
        idx_t nxmax() const { return impl_->nxmax(); }

        /// Number of points per latitude
        const std::vector<idx_t>& nx() const { return impl_->nx(); }

        /// Value of minimum longitude per latitude [default=0]
        const std::vector<double>& xmin() const { return impl_->xmin(); }

        /// Value of maximum longitude per latitude [default=0]
        const std::vector<double>& xmax() const { return impl_->xmax(); }

        /// Value of longitude increment
        const std::vector<double>& dx() const { return impl_->dx(); }

        /// Value of minimum x over entire grid
        double min() const { return impl_->min(); }

        /// Value of maximum x over entire grid
        double max() const { return impl_->max(); }

        Spec spec() const;

        std::string type() const { return impl_->type(); }

    private:
        util::ObjectHandle<XSpace::Implementation> impl_;
    };

    using YSpace = Spacing;

public:
    static std::string static_type();

public:
    Structured(const std::string&, XSpace, YSpace, Projection, Domain);
    Structured(XSpace, YSpace, Projection, Domain);
    Structured(const Structured&, Domain);

    virtual ~Structured() override;

    virtual idx_t size() const override { return npts_; }

    virtual Spec spec() const override;

    /**
   * Human readable name
   * Either the name is the one given at construction as a canonical named grid,
   * or the name "structured"
   */
    virtual std::string name() const override;

    virtual std::string type() const override;

    inline idx_t ny() const { return static_cast<idx_t>(y_.size()); }

    inline idx_t nx(idx_t j) const { return static_cast<idx_t>(nx_[j]); }

    inline idx_t nxmax() const { return nxmax_; }

    inline idx_t nxmin() const { return nxmin_; }

    inline const std::vector<idx_t>& nx() const { return nx_; }

    inline const std::vector<double>& y() const { return y_; }

    inline double dx(idx_t j) const { return dx_[j]; }

    inline double xmin(idx_t j) const { return xmin_[j]; }

    inline double x(idx_t i, idx_t j) const { return xmin_[j] + static_cast<double>(i) * dx_[j]; }

    inline double y(idx_t j) const { return y_[j]; }

    inline void xy(idx_t i, idx_t j, double crd[]) const {
        crd[0] = x(i, j);
        crd[1] = y(j);
    }

    PointXY xy(idx_t i, idx_t j) const { return PointXY(x(i, j), y(j)); }

    PointLonLat lonlat(idx_t i, idx_t j) const { return projection_.lonlat(xy(i, j)); }

    void lonlat(idx_t i, idx_t j, double crd[]) const {
        xy(i, j, crd);
        projection_.xy2lonlat(crd);
    }

    inline bool reduced() const { return nxmax() != nxmin(); }

    bool periodic() const { return periodic_x_; }

    const XSpace& xspace() const { return xspace_; }
    const YSpace& yspace() const { return yspace_; }

    virtual std::unique_ptr<Grid::IteratorXY> xy_begin() const override {
        return std::make_unique<IteratorXY>(*this);
    }
    virtual std::unique_ptr<Grid::IteratorXY> xy_end() const override {
        return std::make_unique<IteratorXY>(*this, false);
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_begin() const override {
        return std::make_unique<IteratorLonLat>(*this);
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_end() const override {
        return std::make_unique<IteratorLonLat>(*this, false);
    }

    gidx_t index(idx_t i, idx_t j) const { return jglooff_[j] + i; }

    void index2ij(gidx_t gidx, idx_t& i, idx_t& j) const {
        if ((gidx < 0) || (gidx >= jglooff_.back())) {
            throw_Exception("Structured::index2ij: gidx out of bounds", Here());
        }
        idx_t ja = 0, jb = jglooff_.size();
        while (jb - ja > 1) {
            idx_t jm = (ja + jb) / 2;
            if (gidx < jglooff_[jm]) {
                jb = jm;
            }
            else {
                ja = jm;
            }
        }
        i = gidx - jglooff_[ja];
        j = ja;
    }

    Config meshgenerator() const override;
    Config partitioner() const override;


protected:  // methods
    virtual void print(std::ostream&) const override;

    virtual void hash(eckit::Hash&) const override;

    virtual RectangularLonLatDomain lonlatBoundingBox() const override;

    void computeTruePeriodicity();

    Domain computeDomain() const;

    void crop(const Domain&);

protected:
    // Minimum number of points across parallels (constant y)
    idx_t nxmin_;

    // Maximum number of points across parallels (constant y)
    idx_t nxmax_;

    /// Total number of unique points in the grid
    idx_t npts_;

    /// Latitude values
    std::vector<double> y_;

    /// Number of points per latitude
    std::vector<idx_t> nx_;

    /// Value of minimum longitude per latitude [default=0]
    std::vector<double> xmin_;

    /// Value of maximum longitude per latitude [default=0]
    std::vector<double> xmax_;

    /// Value of longitude increment
    std::vector<double> dx_;

    /// Periodicity in x-direction
    bool periodic_x_;

    /// Per-row offset
    std::vector<gidx_t> jglooff_;

private:
    std::string name_ = {"structured"};
    XSpace xspace_;
    YSpace yspace_;
    mutable std::string type_;
};

extern "C" {
void atlas__grid__Structured__delete(Structured* This);
const Structured* atlas__grid__Structured(const char* identifier);
const Structured* atlas__grid__Structured__config(util::Config* conf);
const Structured* atlas__grid__regular__RegularGaussian(long N);
Structured* atlas__grid__reduced__ReducedGaussian_int(int nx[], long ny);
Structured* atlas__grid__reduced__ReducedGaussian_long(long nx[], long ny);
Structured* atlas__grid__reduced__ReducedGaussian_int_projection(int nx[], long ny,
                                                                 const Projection::Implementation* projection);
Structured* atlas__grid__reduced__ReducedGaussian_long_projection(long nx[], long ny,
                                                                  const Projection::Implementation* projection);
const Structured* atlas__grid__regular__RegularLonLat(long nx, long ny);
const Structured* atlas__grid__regular__ShiftedLonLat(long nx, long ny);
const Structured* atlas__grid__regular__ShiftedLon(long nx, long ny);
const Structured* atlas__grid__regular__ShiftedLat(long nx, long ny);

void atlas__grid__Structured__nx_array(Structured* This, const idx_t*& nx, idx_t& size);
idx_t atlas__grid__Structured__nx(Structured* This, idx_t j);
idx_t atlas__grid__Structured__ny(Structured* This);
gidx_t atlas__grid__Structured__index(Structured* This, idx_t i, idx_t j);
void atlas__grid__Structured__index2ij(Structured* This, gidx_t gidx, idx_t& i, idx_t& j);
idx_t atlas__grid__Structured__nxmin(Structured* This);
idx_t atlas__grid__Structured__nxmax(Structured* This);
double atlas__grid__Structured__y(Structured* This, idx_t j);
double atlas__grid__Structured__x(Structured* This, idx_t i, idx_t j);
void atlas__grid__Structured__xy(Structured* This, idx_t i, idx_t j, double crd[]);
void atlas__grid__Structured__lonlat(Structured* This, idx_t i, idx_t j, double crd[]);
void atlas__grid__Structured__y_array(Structured* This, const double*& lats, idx_t& size);
int atlas__grid__Structured__reduced(Structured* This);

idx_t atlas__grid__Gaussian__N(Structured* This);
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
