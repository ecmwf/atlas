/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <cstdint>
#include <complex>
#include <string>
#include <map>
#include <memory>

#include "atlas/mdspan.h"

namespace atlas::linalg {

class FFT {
public:
    virtual ~FFT() {}
    virtual const std::string& type() const = 0;

    void plan_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
        do_plan_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
    }
    void plan_inverse_c2r_many(
        mdspan<std::complex<double>,dims<2>,layout_stride> in,
        mdspan<double,dims<2>,layout_stride> out) {
            auto howmany = in.extent(0);
            auto size_in = in.extent(1);
            auto size_out = out.extent(1);
            auto dist_in = in.stride(0);
            auto dist_out = out.stride(0);
        do_plan_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in.data_handle(), out.data_handle());
    }
    void plan_inverse_c2r_many(size_t howmany, size_t size_out, std::complex<double>* in, double* out) const {
        auto size_in = (size_out / 2 ) + 1;
        auto dist_in = size_in;
        auto dist_out = size_out;
        do_plan_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
    }

    void inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
        do_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
    }
    void inverse_c2r_many(
        mdspan<std::complex<double>,dims<2>,layout_stride> in,
        mdspan<double,dims<2>,layout_stride> out) {
            auto howmany = in.extent(0);
            auto size_in = in.extent(1);
            auto size_out = out.extent(1);
            auto dist_in = in.stride(0);
            auto dist_out = out.stride(0);
        do_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in.data_handle(), out.data_handle());
    }
    void inverse_c2r_many(size_t howmany, size_t size_out, std::complex<double>* in, double* out) const {
        auto size_in = (size_out / 2 ) + 1;
        auto dist_in = size_in;
        auto dist_out = size_out;
        inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
    }

    void plan_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) {
        do_plan_inverse_c2r(size_out, in, out);
    }

    void inverse_c2r(size_t size_out, std::complex<double>* in, double* out) {
        do_inverse_c2r(size_out, in, out);
    }

private:

    virtual void do_plan_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const = 0;
    virtual void do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const = 0;
    virtual void do_plan_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const = 0;
    virtual void do_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const = 0;

};

class FFTW : public FFT {
public:
    static const std::string& static_type() {
        static std::string _static_type{"FFTW"};
        return _static_type;
    }
    const std::string& type() const override { return static_type(); }
    FFTW();
    virtual ~FFTW();

private:
    void do_plan_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const override;
    void do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const override;
    void do_plan_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const override;
    void do_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const override;

private:
    struct FFTW_Plans;
    mutable FFTW_Plans* plans_;
};

class pocketfft : public FFT {
public:
    static const std::string& static_type() {
        static std::string _static_type{"pocketfft"};
        return _static_type;
    }
    const std::string& type() const override { return static_type(); }
    virtual ~pocketfft() {}

private:
    void do_plan_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const override;
    void do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const override;
    void do_plan_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const override;
    void do_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const override;
};

} // namespace atlas::linalg
