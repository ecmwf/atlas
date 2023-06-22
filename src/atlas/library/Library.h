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

#include "eckit/filesystem/PathName.h"
#include "eckit/system/Library.h"
#include "eckit/system/Plugin.h"

namespace eckit {
class Parametrisation;
class PathName;
}  // namespace eckit

namespace atlas {

namespace mpi {
void finalise();
void finalize();
}  // namespace mpi

void initialise(int argc, char** argv);
void initialize(int argc, char** argv);
void initialise();
void initialize();
void finalise();
void finalize();

//----------------------------------------------------------------------------------------------------------------------

class Library : public eckit::system::Library {
public:
    static Library& instance();

    virtual std::string version() const override;

    virtual std::string gitsha1(unsigned int count) const override;
    std::string gitsha1() const { return gitsha1(7); }

    void initialise(int argc, char** argv);
    void initialise(const eckit::Parametrisation&);
    void initialise();
    void finalise();

    struct Information {
        friend std::ostream& operator<<(std::ostream& s, const Information& i) {
            i.print(s);
            return s;
        }
        void print(std::ostream&) const;
    };
    Information information() const { return Information(); }

    virtual eckit::Channel& infoChannel() const;
    virtual eckit::Channel& warningChannel() const;
    virtual eckit::Channel& traceChannel() const;
    virtual eckit::Channel& debugChannel() const override;
    bool trace() const { return trace_; }
    virtual bool debug() const override { return debug_; }

    bool traceBarriers() const { return trace_barriers_; }
    bool traceMemory() const { return trace_memory_; }

    Library();

    void registerPlugin(eckit::system::Plugin&);
    void deregisterPlugin(eckit::system::Plugin&);
    const std::vector<eckit::system::Plugin*>& plugins() { return plugins_; }

    std::string cachePath() const;
    std::string dataPath() const;

    std::string linalgDenseBackend() const;
    std::string linalgSparseBackend() const;

    void registerDataPath(const std::string&);

protected:
    virtual const void* addr() const override;

    bool initialized_{false};
    bool debug_{false};
    bool info_{true};
    bool warning_{true};
    bool trace_{false};
    bool trace_memory_{false};
    bool trace_barriers_{false};
    bool trace_report_{false};
    mutable std::unique_ptr<eckit::Channel> info_channel_;
    mutable std::unique_ptr<eckit::Channel> warning_channel_;
    mutable std::unique_ptr<eckit::Channel> trace_channel_;
    mutable std::unique_ptr<eckit::Channel> debug_channel_;

private:
    std::vector<eckit::system::Plugin*> plugins_;
    mutable std::vector<std::string> data_paths_;
    size_t atlas_io_trace_hook_;
};

using Atlas = Library;

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
