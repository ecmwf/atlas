/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <mpi.h>
#include <unordered_map>
#include <vector>

#include "atlas/array/DataType.h"

namespace atlas {
  namespace repartition {
    namespace mpi {

      // Set data type dictionary.
      static const auto mpiDataType =
        std::unordered_map<array::DataType::kind_t, MPI_Datatype>{
          {array::DataType::KIND_REAL64, MPI_DOUBLE},
          {array::DataType::KIND_REAL32, MPI_FLOAT},
          {array::DataType::KIND_INT64, MPI_LONG},
          {array::DataType::KIND_INT32, MPI_INT},
          {array::DataType::KIND_UINT64, MPI_UNSIGNED_LONG},
        };

      /// \brief  Helper class to instantiate and finalise MPI derived types.
      class DataType {

      public:

        /// \brief  Default constructor.
        DataType() = default;

        /// \brief  Constructor.
        DataType(const MPI_Datatype handle);

        /// \brief  Destructor.
        ~DataType();

        /// \brief  Move constructor.
        DataType(DataType&& newDataType);

        /// \brief  Move assignment.
        DataType& operator=(DataType&& newDataType);

        /// \brief  Get data type.
        MPI_Datatype getHandle() const;

        // Disable copying.
        DataType(const DataType&) = delete;
        DataType& operator=(const DataType&) = delete;


      private:

        // MPI data type handle.
        MPI_Datatype handle_{MPI_DATATYPE_NULL};

      };

      /// \brief  MPI Indexed derived type.
      class IndexedType : public DataType {

      public:

        /// \brief  Constructor.
        IndexedType(const std::vector<int>& blockLengths,
          const std::vector<int>& blockDisplacements,
          array::DataType::kind_t arrayDataType);

      };

    }
  }
}
