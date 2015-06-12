/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_Bitflags_h
#define atlas_util_Bitflags_h

#include <string>

namespace atlas {
namespace util {

class Bitflags
{
public:
	static void reset(int& flags, int bit = 0)
	{
		flags = bit;
	}

	static void set(int& flags, int bit)
	{
		flags |= bit;
	}

	static void unset(int& flags, int bit)
	{
		flags &= (~bit);
	}

	static void toggle(int& flags, int bit)
	{
		flags ^= bit;
	}

	static bool check(int flags, int bits)
	{
		return (flags & bits) == bits;
	}

	static bool check_all(int flags, int bits)
	{
		return (flags & bits) == bits;
	}

	static bool check_any(int flags, int bits)
	{
		return flags & bits;
	}

	static std::string bitstr(int flags)
	{
	  char str[9] = {0};
	  int i;
	  for (i=7; i>=0; i--) {
	    str[i] = (flags&1)?'1':'0';
	    flags >>= 1;
	  }
		return std::string(str,9);
	}
};

class Topology : public Bitflags
{
public:
  enum {
    NONE     = 0,
    GHOST    = (1<<1),
    PERIODIC = (1<<2),
    BC       = (1<<3),
    WEST     = (1<<4),
    EAST     = (1<<5),
    NORTH    = (1<<6),
    SOUTH    = (1<<7)
  };
};

} // namespace util
} // namespace atlas

#endif
