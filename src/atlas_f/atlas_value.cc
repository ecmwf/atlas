#include "eckit/value/Value.h"
#include "eckit/value/ListContent.h"
#include "atlas/runtime/ErrorHandling.h"

using eckit::Value;
using eckit::ValueList;
using eckit::ListContent;
using eckit::makeVectorValue;

namespace atlas {

extern "C" {
Value* atlas__Value__new_int    (int         value) { return new Value(value); }
Value* atlas__Value__new_long   (long        value) { return new Value(value); }
Value* atlas__Value__new_float  (float       value) { return new Value(value); }
Value* atlas__Value__new_double (double      value) { return new Value(value); }
Value* atlas__Value__new_string (const char* value) { return new Value(value); }

Value* atlas__Value__new_array_int (int value[], int size)
{
  std::vector<int> v(value,value+size);
  return new Value( makeVectorValue(v) );
}

Value* atlas__Value__new_array_long (long value[], int size)
{
  std::vector<long> v(value,value+size);
  return new Value( makeVectorValue(v) );
}

Value* atlas__Value__new_array_float (float value[], int size)
{
  std::vector<float> v(value,value+size);
  return new Value( makeVectorValue(v) );
}

Value* atlas__Value__new_array_double (double value[], int size)
{
  std::vector<double> v(value,value+size);
  return new Value( makeVectorValue(v) );
}

void atlas__Value__delete (Value* This) { delete This; }

void atlas__Value__int (Value* This, int &value) { value = *This; }
void atlas__Value__long (Value* This, long &value) { value = *This; }
void atlas__Value__float (Value* This, float &value) { value = double(*This); }
void atlas__Value__double (Value* This, double &value) { value = *This; }

void atlas__Value__string (Value* This, char* &value, int &size, int &allocated)
{
  std::string tmp = *This;
  size = tmp.size();
  value = new char[tmp.size()+1];
  strcpy(value,tmp.c_str());
  allocated = true;
}

void atlas__Value__array_int (Value* This, int* &value, int &size, int &allocated)
{
  std::vector<eckit::Value> tmp = *This;
  size = tmp.size();
  value = new int[size];
  for(int j = 0; j < size; ++j) value[j] = tmp[j];
  allocated = true;
}

void atlas__Value__array_long (Value* This, long* &value, int &size, int &allocated)
{
  std::vector<eckit::Value> tmp = *This;
  size = tmp.size();
  value = new long[size];
  for(int j = 0; j < size; ++j) value[j] = tmp[j];
  allocated = true;
}

void atlas__Value__array_float (Value* This, float* &value, int &size, int &allocated)
{
  std::vector<eckit::Value> tmp = *This;
  size = tmp.size();
  value = new float[size];
  for(int j = 0; j < size; ++j) value[j] = double(tmp[j]);
  allocated = true;
}

void atlas__Value__array_double (Value* This, double* &value, int &size, int &allocated)
{
  std::vector<eckit::Value> tmp = *This;
  size = tmp.size();
  value = new double[size];
  for(int j = 0; j < size; ++j) value[j] = tmp[j];
  allocated = true;
}

}

}
