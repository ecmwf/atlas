/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_testing_UnitTest_hpp
#define atlas_testing_UnitTest_hpp

/// @file UnitTest.hpp
/// @author Willem Deconinck
/// This file implements a simple UnitTest class
/// Example usage:
/**
  class MyUnitTest : public UnitTest
  {
  public:
    MyUnitTest(int argc, char *argv[]) : UnitTest(argc,argv) {}
    virtual void run_tests()
    { 
      // Implementation here
      ATLAS_CHECK_EQUAL( 1, 2 );
    }
  };

  int main(int argc, char *argv[])
  {
    MyUnitTest tests(argc,argv);
    tests.start();
    return tests.exit_status();
  } 
*/

#include "atlas/testing/Macros.hpp"
#include "atlas/mpl/MPL.hpp"

namespace atlas {
class UnitTest
{
public:
  UnitTest(int argc, char *argv[]) : argc_(argc), argv_(argv)
  {
    failed_tests_=0;
    MPL::init(argc_,argv_);
  }
  
  ~UnitTest()
  {
    MPL::finalize();
  }
  
  void start() {
    try{
      run_tests();
    }
    catch( std::exception& e )
    {
      std::cout << "Exiting with exception\n" << e.what() << std::endl;
    }
  }
  
  int exit_status() { return failed_tests_ > 0 ? 1 : 0; }

protected:
    
  virtual void run_tests()=0;
    
  template< typename T1, typename T2 >
  int assert_equal(const T1& v1, const T2& v2)
  {
    if( v1 != v2 ){
      std::stringstream ss; ss << "[" << v1 << " != " << v2 << "]" << std::endl;
      error_msg_ = ss.str();
      ++failed_tests_;
      return ATLAS_CHECK_FAILED;
    }
    return ATLAS_CHECK_SUCCESS;
  }
  
  int assert( bool assertion )
  {
    if( ! assertion ) {
      error_msg_.clear();
      ++failed_tests_;
      return ATLAS_CHECK_FAILED;
    }
    return ATLAS_CHECK_SUCCESS;
  }

protected:
  int argc_;
  char **argv_;
  std::string error_msg_;
  int failed_tests_;
};

}

#endif
