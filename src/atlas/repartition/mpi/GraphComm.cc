/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/repartition/mpi/GraphComm.h"

namespace atlas {
  namespace repartition {
    namespace mpi {

      // Constructor.
      GraphComm::GraphComm(const atlas::mpi::Comm& oldComm,
        const std::vector<int>& sources, const std::vector<int>& destinations) {

        const auto commHandle = oldComm.communicator();
        const auto inDegree = static_cast<int>(sources.size());
        const auto outDegree = static_cast<int>(destinations.size());

        CALL_MPI(MPI_Dist_graph_create_adjacent(
          commHandle, inDegree, sources.data(), MPI_UNWEIGHTED, outDegree,
          destinations.data(), MPI_UNWEIGHTED, MPI_INFO_NULL, false, &handle_));

        return;
      }

      // Destructor.
      GraphComm::~GraphComm() {
        if (handle_ != MPI_COMM_NULL) CALL_MPI(MPI_Comm_free(&handle_));
        return;
      }

      // Move constructor.
      GraphComm::GraphComm(GraphComm&& newGraphComm) {

        std::swap(handle_, newGraphComm.handle_);

        return;
      }

      // Move assignment operator.
      GraphComm& GraphComm::operator=(GraphComm &&newGraphComm) {

        std::swap(handle_, newGraphComm.handle_);

        return *this;
      }

      void GraphComm::neighbourAllToAllw(const void* const sendBuffer,
        const std::vector<int>& sendCounts,
        const std::vector<nByte>& sendDisplacements,
        const std::vector<DataType>& sendDataTypes,
        void* const recvBuffer,
        const std::vector<int>& recvCounts,
        const std::vector<nByte>& recvDisplacements,
        const std::vector<DataType>& recvDataTypes) const {

        // Get data type handles.
        auto getDataTypeHandles = [](const std::vector<DataType>& dataTypes){
          auto dataTypeHandles = std::vector<MPI_Datatype>(dataTypes.size());
          std::transform(dataTypes.cbegin(), dataTypes.cend(),
          dataTypeHandles.begin(), [](const DataType& elem){
            return elem.getHandle();
          });
          return dataTypeHandles;
        };

        auto sendDataTypeHandles = getDataTypeHandles(sendDataTypes);
        auto recvDataTypeHandles = getDataTypeHandles(recvDataTypes);

        CALL_MPI(MPI_Neighbor_alltoallw(sendBuffer, sendCounts.data(),
          sendDisplacements.data(), sendDataTypeHandles.data(), recvBuffer,
          recvCounts.data(), recvDisplacements.data(),
          recvDataTypeHandles.data(), handle_));

        return;
      }
    }
  }
}
