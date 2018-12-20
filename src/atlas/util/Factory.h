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
#include <vector>

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace util {

class FactoryBase;

class FactoryRegistry {
protected:
    FactoryRegistry( const std::string& factory );
    ~FactoryRegistry();

private:
    mutable std::mutex mutex_;
    std::map<std::string, FactoryBase*> factories_;
    std::string factory_;

public:
    std::vector<std::string> keys() const;
    void list( std::ostream& ) const;
    bool has( const std::string& builder ) const;
    void add( const std::string& builder, FactoryBase* );
    void remove( const std::string& builder );
    FactoryBase* get( const std::string& builder ) const;
};

template <typename T>
struct FactoryRegistryT : public FactoryRegistry {
public:
    static FactoryRegistryT<T>& instance() {
        static FactoryRegistryT<T> env( T::className() );
        return env;
    }

private:
    FactoryRegistryT( const std::string& factory ) : FactoryRegistry( factory ) {}
};

class FactoryBase {
private:
    FactoryRegistry& registry_;
    std::string builder_;

protected:
    FactoryBase( FactoryRegistry&, const std::string& builder );
    virtual ~FactoryBase();
    friend class FactoryRegistry;
};

template <typename T>
class Factory : public FactoryBase {
public:
    static std::vector<std::string> keys() { return registry().keys(); }
    static void list( std::ostream& out ) { return registry().list( out ); }
    static bool has( const std::string& builder ) { return registry().has( builder ); }
    static T* get( const std::string& builder ) { return dynamic_cast<T*>( registry().get( builder ) ); }

    Factory( const std::string& builder = "" ) : FactoryBase( registry(), builder ) {}

protected:
    virtual ~Factory() = default;
    static FactoryRegistryT<T>& registry() { return FactoryRegistryT<T>::instance(); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
