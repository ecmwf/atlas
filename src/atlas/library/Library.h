/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>
#include <string>
#include "eckit/system/Library.h"

namespace eckit {
  class Parametrisation;
}

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class Library : public eckit::system::Library {
public:

    Library();

    static Library& instance();

    virtual std::string version() const override;

    virtual std::string gitsha1(unsigned int count) const override;
    std::string gitsha1() const { return gitsha1(7); }

    void initialise(int argc, char **argv);
    void initialise(const eckit::Parametrisation&);
    void initialise();
    void finalise();
    
    struct Info {
      friend std::ostream& operator<<(std::ostream& s, const Info& i) { i.print(s); return s; }
      void print(std::ostream&) const;
    };
    Info info() const { return Info(); }

    class Timer {
    public:
        bool barriers() const;
        std::ostream& channel() const;
    private:
        friend class Library;
        bool barriers_{false};
        std::ostream* channel_;
    };

    const Timer& timer() const { return timer_; }

protected:

    virtual const void* addr() const override;

    Timer timer_;

};

typedef Library Atlas;

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

