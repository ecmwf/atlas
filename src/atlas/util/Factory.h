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
#include <memory>
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
    FactoryRegistry(const std::string& factory);
    virtual ~FactoryRegistry();

private:
    mutable std::mutex mutex_;
    std::string factory_;
    std::map<std::string, FactoryBase*> factories_;

public:
    const std::string& factory() const { return factory_; }
    std::vector<std::string> keys() const;
    void list(std::ostream&) const;
    bool has(const std::string& builder) const;
    void add(const std::string& builder, FactoryBase*);
    void remove(const std::string& builder);
    FactoryBase* get(const std::string& builder) const;
};

template <typename T>
struct FactoryRegistryT : public FactoryRegistry {
public:
    static std::shared_ptr<FactoryRegistryT<T>> instance() {
        static std::shared_ptr<FactoryRegistryT<T>> env(new FactoryRegistryT<T>(T::className()));
        return env;
    }
    virtual ~FactoryRegistryT() {}

protected:
    FactoryRegistryT(const std::string& factory): FactoryRegistry(factory) {}
};

class FactoryBase {
private:
    FactoryRegistry& registry_;
    std::string builder_;
    std::shared_ptr<FactoryRegistry> attached_registry_;

protected:
    FactoryBase(FactoryRegistry&, const std::string& builder);
    virtual ~FactoryBase();
    void attach_registry(const std::shared_ptr<FactoryRegistry>& registry) { attached_registry_ = registry; }
    friend class FactoryRegistry;

public:
    const std::string& factoryBuilder() const { return builder_; }
    const std::string& factoryName() const { return registry_.factory(); }
};

template <typename T>
class Factory : public FactoryBase {
public:
    static std::vector<std::string> keys() { return registry().keys(); }
    static void list(std::ostream& out) { return registry().list(out); }
    static bool has(const std::string& builder) { return registry().has(builder); }
    static T* get(const std::string& builder) { return dynamic_cast<T*>(registry().get(builder)); }

    Factory(const std::string& builder = ""): FactoryBase(registry(), builder) {
        if (not builder.empty()) {
            attach_registry(FactoryRegistryT<T>::instance());
        }
    }

protected:
    virtual ~Factory(){};
    static FactoryRegistry& registry() { return *FactoryRegistryT<T>::instance().get(); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
