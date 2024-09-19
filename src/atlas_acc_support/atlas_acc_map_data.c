/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


// This file needs to be compiled with an OpenACC capable C compiler that is
// compatible with the Fortran compiler

#ifndef _OPENACC
#error atlas_acc_map_data must be compiled with OpenACC capable compiler
#endif

#include <openacc.h>
#include <string.h>

#include "atlas_acc_map_data.h"

void atlas_acc_map_data(void* cpu_ptr, void* gpu_ptr, unsigned long bytes) {
  acc_map_data(cpu_ptr, gpu_ptr, bytes);
}


void atlas_acc_unmap_data(void* cpu_ptr) {
  acc_unmap_data(cpu_ptr);
}

int atlas_acc_is_present(void* cpu_ptr, unsigned long bytes) {
  return acc_is_present(cpu_ptr, bytes);
}

void* atlas_acc_deviceptr(void* cpu_ptr) {
  return acc_deviceptr(cpu_ptr);
}

atlas_acc_device_t atlas_acc_get_device_type() {
  acc_device_t device_type = acc_get_device_type();
  if (device_type == acc_device_host || device_type == acc_device_none) {
    return atlas_acc_device_host;
  }
  return atlas_acc_device_not_host;
}

#define BUFFER_PRINTF(buffer, ...) snprintf(buffer + strlen(buffer), sizeof(buffer) - strlen(buffer), __VA_ARGS__);

const char* atlas_acc_version_str() {
  static char buffer[16];
  if( strlen(buffer) != 0 ) {
    return buffer;
  }

  BUFFER_PRINTF(buffer, "%i", _OPENACC);
  switch( _OPENACC ) {
    case 201111: BUFFER_PRINTF(buffer, " (1.0)"); break;
    case 201306: BUFFER_PRINTF(buffer, " (2.0)"); break;
    case 201308: BUFFER_PRINTF(buffer, " (2.0)"); break;
    case 201510: BUFFER_PRINTF(buffer, " (2.5)"); break;
    case 201711: BUFFER_PRINTF(buffer, " (2.6)"); break;
    default: break;
  }
  return buffer;
}
const char* atlas_acc_info_str() {
    static char buffer[1024];
    if( strlen(buffer) != 0 ) {
      return buffer;
    }
    // Platform information
    acc_device_t devicetype = acc_get_device_type();
    int num_devices = acc_get_num_devices(devicetype);
    BUFFER_PRINTF(buffer, "    OpenACC version: %s\n", atlas_acc_version_str());
    if (devicetype == acc_device_host ) {
      BUFFER_PRINTF(buffer, "    No OpenACC GPU devices available\n");
      return buffer;
    }
    BUFFER_PRINTF(buffer, "    Number of OpenACC devices: %i\n", num_devices);
    int device_num  = acc_get_device_num(devicetype);
    BUFFER_PRINTF(buffer, "    OpenACC Device number: %i\n", device_num);
    BUFFER_PRINTF(buffer, "    OpenACC acc_device_t (enum value): %d", devicetype);
    switch( devicetype ) {
      case acc_device_host: BUFFER_PRINTF(buffer, " (acc_device_host)"); break;
      case acc_device_none: BUFFER_PRINTF(buffer, " (acc_device_none)"); break;
      case acc_device_not_host: BUFFER_PRINTF(buffer, " (acc_device_not_host)"); break;
      case acc_device_nvidia: BUFFER_PRINTF(buffer, " (acc_device_nvidia)"); break;
      default: break;
    }
    BUFFER_PRINTF(buffer, "\n");

    // acc_get_property, acc_get_property_string introduced with OpenACC 2.6
#if _OPENACC >= 201711
    long int    mem             = acc_get_property(       device_num, devicetype, acc_property_memory);
    long int    free_mem        = acc_get_property(       device_num, devicetype, acc_property_free_memory);
    const char *property_name   = acc_get_property_string(device_num, devicetype, acc_property_name);
    const char *property_vendor = acc_get_property_string(device_num, devicetype, acc_property_vendor );
    const char *property_driver = acc_get_property_string(device_num, devicetype, acc_property_driver );
    // if (mem > 0)
    {
      BUFFER_PRINTF(buffer, "    Memory on OpenACC device: %li\n", mem);
    }
    // if (free_mem > 0)
    {
       BUFFER_PRINTF(buffer, "    Free Memory on OpenACC device: %li\n", free_mem);
    }
    // if (property_name != NULL)
    {
      BUFFER_PRINTF(buffer, "    OpenACC device name: %s\n", property_name);
    }
    // if (property_vendor != NULL)
    {
      BUFFER_PRINTF(buffer, "    OpenACC device vendor: %s\n", property_vendor);
    }
    // if (property_driver != NULL)
    {
      BUFFER_PRINTF(buffer, "    OpenACC device driver: %s\n", property_driver);
    }
#endif
    return buffer;
}

int atlas_acc_get_num_devices() {
  return acc_get_num_devices(acc_get_device_type());
}
