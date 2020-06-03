
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
               const atlas::Field & f, size_t ngptot) 
               : _v (v), _ind (ind), _f (f), _ngptot (ngptot) 
  {
    _ngpblks = _v.shape (0);
    _nproma  = _v.shape (1);
    _dlength = _v.shape (2);
    if (_ngptot == 0)
      _ngptot = _ngpblks * _nproma;
    _size = _ngptot * _dlength;
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

  size_t ngpblks () const
  {
    return _ngpblks;
  }

  size_t ngptot () const
  {
    return _ngptot;
  }

  size_t nproma () const
  {
    return _nproma;
  }

  size_t dlength () const
  {
    return _dlength;
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
  size_t _ngptot;   // Total number of elements in a distributed field
  size_t _nproma;   // Blocking factor; may be equal to NGPTOT
  size_t _ngpblks;  // Number of blocks
  size_t _size;     // Size in bytes
  size_t _dlength;  // Element length
};

void createIoFieldDescriptorsBlocked 
  (atlas::Field & f, std::vector<ioFieldDesc> & list, atlas::idx_t bdim, atlas::idx_t gdim = -1, size_t ngptot = 0);

void createIoFieldDescriptors (atlas::Field & f, std::vector<ioFieldDesc> & list, size_t ngptot = 0, atlas::idx_t gdim = -1);
void createIoFieldDescriptors (atlas::FieldSet & s, std::vector<ioFieldDesc> & list, size_t ngptot = 0, atlas::idx_t gdim = -1);
