// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "NumericalAlgorithms/Convergence/HasConverged.hpp"
#include "Utilities/PrettyType.hpp"

namespace Parallel::Tags {

/*!
 * \brief Holds an `IterationId` that identifies a step in a parallel algorithm
 */
template <typename Label>
struct IterationId : db::SimpleTag {
  static std::string name() noexcept {
    return "IterationId(" + pretty_type::short_name<Label>() + ")";
  }
  using type = size_t;
};

/*!
 * \brief Holds a `Convergence::HasConverged` flag that signals the parallel
 * algorithm has converged, along with the reason for convergence.
 */
template <typename Label>
struct HasConverged : db::SimpleTag {
  static std::string name() noexcept {
    return "HasConverged(" + pretty_type::short_name<Label>() + ")";
  }
  using type = Convergence::HasConverged;
};

}  // namespace Parallel::Tags
