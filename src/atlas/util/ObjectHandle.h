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

namespace atlas {
namespace util {

class Object;

class ObjectHandleBase {
public:
    ObjectHandleBase() = default;
    ObjectHandleBase( const Object* object ) : object_( const_cast<Object*>( object ) ) { attach(); }

    virtual ~ObjectHandleBase() {
        if ( object_ ) {
            release();
        }
    }

    const ObjectHandleBase& operator=( const ObjectHandleBase& other ) {
        if ( object_ != other.object_ ) {
            assign( other );
        }
        return *this;
    }

    operator bool() const { return object_ != nullptr; }

    void reset( const Object* other ) {
        if ( object_ != other ) {
            assign( other );
        }
    }

    int owners() const;

private:
    void release();

    void assign( const ObjectHandleBase& other );

    void assign( const Object* other );

    void attach();

    bool null() const { return ( object_ == nullptr ); }

protected:
    Object* object_{nullptr};
};

template <typename T>
class ObjectHandle : public ObjectHandleBase {
public:
    using Implementation = T;
    using Handle         = ObjectHandle<T>;

public:
    T* get() { return reinterpret_cast<T*>( object_ ); }
    const T* get() const { return reinterpret_cast<const T*>( object_ ); }
    ObjectHandle() = default;
    ObjectHandle( const T* object ) : ObjectHandleBase( reinterpret_cast<const Object*>( object ) ) {}
    ObjectHandle( const ObjectHandle& handle ) : ObjectHandleBase( reinterpret_cast<const Object*>( handle.get() ) ) {}
    const T* operator->() const { return get(); }
    T* operator->() { return get(); }
    const T& operator*() const { return *get(); }
    T& operator*() { return *get(); }
    void reset( const T* object ) { ObjectHandleBase::reset( reinterpret_cast<const Object*>( object ) ); }
};

}  // namespace util
}  // namespace atlas
