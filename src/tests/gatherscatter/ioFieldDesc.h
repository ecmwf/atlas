
#pragma once

#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>


typedef char byte;

class ioFieldDesc
{
public:
  ioFieldDesc (atlas::array::ArrayView<byte,3> & v, 
               const std::vector<atlas::idx_t> & ind,
               const atlas::Field & f, size_t ldim) 
               : _v (v), _ind (ind), _f (f), _ldim (ldim) 
  {
    _nblk = _v.shape (0);
    _dlen = _v.shape (2);
    if (_ldim == 0)
      _ldim = _v.shape (1);
    _size = _nblk * _ldim * _dlen;
  }


  atlas::array::ArrayView<byte,3> & view ()
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

  size_t nblk () const
  {
    return _nblk;
  }

  size_t ldim () const
  {
    return _ldim;
  }

  size_t dlen () const
  {
    return _dlen;
  }

  byte & operator () (int jblk, int jlon, int k)
  {
    return _v (jblk, jlon, k);
  }

  byte operator () (int jblk, int jlon, int k) const
  {
    return _v (jblk, jlon, k);
  }

private:
  atlas::array::ArrayView<byte,3> _v; // NGPBKLS, NPROMA, sizeof (element)
  const std::vector<atlas::idx_t> _ind;
  const atlas::Field & _f;
  int _owner = 0;
  size_t _ldim;
  size_t _nblk;
  size_t _size;
  size_t _dlen;
};

void createIoFieldDescriptors (atlas::Field & f, std::vector<ioFieldDesc> & list, size_t ldim = 0, atlas::idx_t gdim = -1);
void createIoFieldDescriptors (atlas::FieldSet & s, std::vector<ioFieldDesc> & list, size_t ldim = 0, atlas::idx_t gdim = -1);
