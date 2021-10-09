/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/redistribution/detail/RedistributionImpl.h"

namespace atlas {
namespace redistribution {
namespace detail {

class RedistributeGeneric : public RedistributionImpl {
public:

  static std::string static_type() { return "RedistributeGeneric"; }

  std::string type() const override { return static_type(); }

  void setup( const FunctionSpace& sourceFunctionSpace,
              const FunctionSpace& targetFunctionSpace) override;

  void execute( const Field& sourceField, Field& targetField ) const override;

  void execute( const FieldSet& sourceFieldSet, FieldSet& targetFieldSet ) const override;

private:

  // Determine datatype.
  void doExecute( const Field& sourceField, Field& targetField ) const;
  // Determine rank.
  template <typename Value>
  void doExecute( const Field& sourceField, Field& targetField ) const;
  // Perform redistribution.
  template <typename Value, int Rank>
  void doExecute( const Field& sourceField, Field& targetField ) const;

  // Local indices to send to each PE
  std::vector<idx_t> sourceLocalIdx_{};

  // Local indices to receive from each PE.
  std::vector<idx_t> targetLocalIdx_{};

  // Partial sum of number of columns to send to each PE.
  std::vector<int> sourceDisps_{};

  // Partial sum of number of columns to receive from each PE.
  std::vector<int> targetDisps_{};

};

} // namespace detail
} // namespace redistribution
} // namespace atlas
