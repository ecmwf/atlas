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

#include <iosfwd>
#include <memory>
#include <string>

#include "eckit/system/Library.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {

namespace mpi {
void finalise();
void finalize();
}  // namespace mpi

void initialise( int argc, char** argv );
void initialize( int argc, char** argv );
void initialise();
void initialize();
void finalise();
void finalize();

//----------------------------------------------------------------------------------------------------------------------

class Library : public eckit::system::Library {
public:
    static Library& instance();

    virtual std::string version() const override;

    virtual std::string gitsha1( unsigned int count ) const override;
    std::string gitsha1() const { return gitsha1( 7 ); }

    void initialise( int argc, char** argv );
    void initialise( const eckit::Parametrisation& );
    void initialise();
    void finalise();

    struct Information {
        friend std::ostream& operator<<( std::ostream& s, const Information& i ) {
            i.print( s );
            return s;
        }
        void print( std::ostream& ) const;
    };
    Information information() const { return Information(); }

    virtual eckit::Channel& infoChannel() const;
    virtual eckit::Channel& warningChannel() const;
    virtual eckit::Channel& traceChannel() const;
    virtual eckit::Channel& debugChannel() const override;
    bool trace() const { return trace_; }
    virtual bool debug() const override { return debug_; }

    bool traceBarriers() const { return trace_barriers_; }

    Library();

protected:
    virtual const void* addr() const override;

    bool debug_{false};
    bool info_{true};
    bool warning_{true};
    bool trace_{false};
    bool trace_barriers_{false};
    bool trace_report_{false};
    mutable std::unique_ptr<eckit::Channel> info_channel_;
    mutable std::unique_ptr<eckit::Channel> warning_channel_;
    mutable std::unique_ptr<eckit::Channel> trace_channel_;
    mutable std::unique_ptr<eckit::Channel> debug_channel_;
};

using Atlas = Library;

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
