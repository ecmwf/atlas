/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/repartition/mpi/call.h"
#include "atlas/repartition/mpi/DataType.h"

namespace atlas {
  namespace repartition {
    namespace mpi {

      // Data type constructor.
      DataType::DataType(MPI_Datatype handle) : handle_(handle) {
        CALL_MPI(MPI_Type_commit(&handle));
        return;
      }

      // Data type destructor.
      DataType::~DataType() {
        if(handle_ != MPI_DATATYPE_NULL) CALL_MPI(MPI_Type_free(&handle_));
        return;
      }

      // Move constructor.
      DataType::DataType(DataType&& newDataType) {

        std::swap(handle_, newDataType.handle_);

        return;
      }

      // Move assignment operator.
      DataType& DataType::operator=(DataType &&newDataType) {

        std::swap(handle_, newDataType.handle_);

        return *this;
      }


      // Get data type handle.
      MPI_Datatype DataType::getHandle() const {
        return handle_;
      }

      // Indexed type constructor.
      IndexedType::IndexedType(const std::vector<int>& blockLengths,
        const std::vector<int>& blockDisplacements,
        array::DataType::kind_t arrayDataType) : DataType(
          [&](){
            // Declare result.
            MPI_Datatype newDataType{};

            // Create indexed data type.
            auto oldDataType = mpiDataType.at(arrayDataType);
            const auto count = static_cast<int>(blockLengths.size());
            CALL_MPI(MPI_Type_indexed(count, blockLengths.data(),
              blockDisplacements.data(), oldDataType, &newDataType));

            return newDataType;
        }()) {}

    }
  }
}
