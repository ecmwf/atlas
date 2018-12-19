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

#include <map>
#include <mutex>
#include <string>

namespace atlas {
namespace util {

class ObjectFactory;

class ObjectFactoryRegistry {
protected:
    ObjectFactoryRegistry( const std::string& factory );
    ~ObjectFactoryRegistry();

private:
    mutable std::mutex mutex_;
    std::map<std::string, ObjectFactory*> factories_;
    std::string factory_;

public:
    void list( std::ostream& ) const;
    bool has( const std::string& builder ) const;
    void add( const std::string& builder, ObjectFactory* );
    void remove( const std::string& builder );
    ObjectFactory* get( const std::string& builder ) const;
};

template <typename T>
struct ObjectFactoryRegistryT : public ObjectFactoryRegistry {
public:
    static ObjectFactoryRegistryT<T>& instance() {
        static ObjectFactoryRegistryT<T> env( T::classname() );
        return env;
    }

private:
    ObjectFactoryRegistryT( const std::string& factory ) : ObjectFactoryRegistry( factory ) {}
};

class ObjectFactory {
private:
    ObjectFactoryRegistry& registry_;
    std::string builder_;

protected:
    ObjectFactory( ObjectFactoryRegistry&, const std::string& builder );
    virtual ~ObjectFactory();
    friend class ObjectFactoryRegistry;
};

template <typename T>
class ObjectFactoryT : public ObjectFactory {
public:
    using Factory = ObjectFactoryT<T>;

public:
    static void list( std::ostream& out ) { return registry().list( out ); }
    static bool has( const std::string& builder ) { return registry().has( builder ); }
    static T* get( const std::string& builder ) { return dynamic_cast<T*>( registry().get( builder ) ); }

    ObjectFactoryT( const std::string& builder = "" ) : ObjectFactory( registry(), builder ) {}

protected:
    virtual ~ObjectFactoryT() = default;
    static ObjectFactoryRegistryT<T>& registry() { return ObjectFactoryRegistryT<T>::instance(); }
};


}  // namespace util
}  // namespace atlas
