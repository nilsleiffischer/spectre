// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/FixedHashMap.hpp"
#include "Domain/Direction.hpp"
#include "Domain/ElementId.hpp"
#include "Domain/MaxNumberOfNeighbors.hpp"
#include "Options/Options.hpp"
#include "Utilities/TMPL.hpp"

namespace LinearSolver {
namespace schwarz_detail {

namespace OptionTags {

template <typename OptionsGroup>
struct Overlap {
  using type = size_t;
  using group = OptionsGroup;
  static constexpr OptionString help =
      "Number of points a subdomain overlaps with its neighbor";
};

template <typename SolverType, typename OptionsGroup>
struct SubdomainSolver {
  using type = SolverType;
  using group = OptionsGroup;
  static constexpr OptionString help =
      "Options for the linear solver on subdomains";
};

}  // namespace OptionTags

namespace Tags {

template <typename OptionsGroup>
struct Overlap : db::SimpleTag {
  static std::string name() noexcept {
    return option_name<OptionsGroup>() + "(Overlap)";
  }
  using type = size_t;

  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::Overlap<OptionsGroup>>;
  static type create_from_options(const type& value) noexcept {
    return value;
  }
};

template <typename OptionsGroup>
struct SubdomainSolverBase : db::BaseTag {};

template <typename SolverType, typename OptionsGroup>
struct SubdomainSolver : SubdomainSolverBase<OptionsGroup>, db::SimpleTag {
  static std::string name() noexcept {
    return option_name<OptionsGroup>() + "(SubdomainSolver)";
  }
  using type = SolverType;

  static constexpr bool pass_metavariables = false;
  using option_tags =
      tmpl::list<OptionTags::SubdomainSolver<SolverType, OptionsGroup>>;
  static type create_from_options(const type& value) noexcept {
    return value;
  }
};

template <typename SubdomainOperator, typename OptionsGroup>
struct SubdomainBoundaryData : db::SimpleTag {
  static std::string name() noexcept {
    return option_name<OptionsGroup>() + "(SubdomainBoundaryData)";
  }
  using type = typename SubdomainOperator::SubdomainDataType::BoundaryDataType;
};

}  // namespace Tags
}  // namespace schwarz_detail
}  // namespace LinearSolver
