#pragma once

#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class LonLatProjectionT : public ProjectionImpl {
private:
    friend class ProjectionImpl;
    LonLatProjectionT() : ProjectionImpl() {}

public:
    // constructor
    LonLatProjectionT( const eckit::Parametrisation& );

    // destructor
    ~LonLatProjectionT() {}

    // class name
    static std::string static_type() { return Rotation::typePrefix() + "lonlat"; }
    virtual std::string type() const override { return static_type(); }

    // projection and inverse projection
    virtual void xy2lonlat( double crd[] ) const override { rotation_.rotate( crd ); }
    virtual void lonlat2xy( double crd[] ) const override { rotation_.unrotate( crd ); }

    virtual bool strictlyRegional() const override { return false; }

    // specification
    virtual Spec spec() const override;

    virtual std::string units() const override { return "degrees"; }

    virtual operator bool() const override { return rotation_.rotated(); }

    virtual void hash( eckit::Hash& ) const override;

private:
    Rotation rotation_;
};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated> RotatedLonLatProjection;

// --------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
