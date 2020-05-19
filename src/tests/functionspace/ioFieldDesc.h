
#pragma once

#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>


typedef char byte;

class ioFieldDesc
{
public:
  ioFieldDesc (atlas::array::ArrayView<byte,2> & v, 
               const std::vector<atlas::idx_t> & ind,
               const atlas::Field & f, size_t ldim) 
               : _v (v), _ind (ind), _f (f), _ldim (ldim) 
  {
    _dlen = _v.shape (1);
    if (_ldim == 0)
      _ldim = _v.shape (0);
    _size = _ldim * _v.shape (1);
  }

  atlas::array::ArrayView<byte,2> & view ()
  {
    return _v;
  }

  const std::vector<atlas::idx_t> & indices () const
  {
    return _ind;
  }

  const atlas::Field & field () const
  {
    return _f;
  }

  int & owner () 
  {
    return _owner;
  }

  int owner () const
  {
    return _owner;
  }

  size_t size () const
  {
    return _size;
  }

  size_t dlen () const
  {
    return _dlen;
  }

  void pack (byte * buffer) const
  {
    for (int i = 0; i < _ldim; i++)
    for (int j = 0; j < _dlen; j++)
      buffer[i*_dlen+j] = _v (i, j);
  }
 
  byte & operator () (int i, int j)
  {
    return _v (i, j);
  }

private:
  atlas::array::ArrayView<byte,2> _v;
  const std::vector<atlas::idx_t> _ind;
  const atlas::Field & _f;
  int _owner = 0;
  size_t _ldim;
  size_t _size;
  size_t _dlen;
};

void createIoFieldDescriptors (atlas::Field & f, std::vector<ioFieldDesc> & list, size_t ldim = 0);
void createIoFieldDescriptors (atlas::FieldSet & s, std::vector<ioFieldDesc> & list, size_t ldim = 0);
