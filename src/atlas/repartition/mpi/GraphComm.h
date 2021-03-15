/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/repartition/mpi/DataType.h"

namespace atlas {
  namespace repartition {
    namespace mpi {

      // MPI byte displacement type.
      using nByte = MPI_Aint;

      /// \brief  MPI graph communicator class.
      class GraphComm {

      public:

        /// \brief  Default onstructor.
        GraphComm() = default;

        /// \brief  Constructor.
        GraphComm(const atlas::mpi::Comm& oldComm,
          const std::vector<int>& sources,
          const std::vector<int>& destinations);

        /// \brief  Destructor.
        ~GraphComm();

        /// \brief  Move constructor.
        GraphComm(GraphComm&& newGraphComm);

        /// \brief  Move assignment.
        GraphComm& operator=(GraphComm&& newGraphComm);

        /// \brief  Neighbour all-to-all with variable data types.
        void neighbourAllToAllw(const void* const sendBuffer,
          const std::vector<int>& sendCounts,
          const std::vector<nByte>& sendDisplacements,
          const std::vector<DataType>& sendDataTypes,
          void* const recvBuffer,
          const std::vector<int>& recvCounts,
          const std::vector<nByte>& recvDisplacements,
          const std::vector<DataType>& recvDataTypes) const;

        // Disable copying.
        GraphComm(const GraphComm&) = delete;
        GraphComm& operator=(const GraphComm&) = delete;

      private:

        // Communicator handle.
        MPI_Comm handle_{MPI_COMM_NULL};

      };

    }
  }
}
