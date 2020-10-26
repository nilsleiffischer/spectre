// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Subitems.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Options/Options.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/OverlapHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace LinearSolver::Schwarz {

/// Option tags related to the Schwarz solver
namespace OptionTags {

template <typename OptionsGroup>
struct MaxOverlap {
  using type = size_t;
  using group = OptionsGroup;
  static constexpr Options::String help =
      "Number of points that subdomains can extend into neighbors";
};

template <typename SolverType, typename OptionsGroup>
struct SubdomainSolver {
  using type = SolverType;
  using group = OptionsGroup;
  static constexpr Options::String help =
      "Options for the linear solver on subdomains";
};

template <typename OptionsGroup>
struct SubdomainPreconditionerGroup {
  static std::string name() noexcept { return "SubdomainPreconditioner"; }
  using group = OptionsGroup;
  static constexpr Options::String help =
      "The preconditioner for subdomain solves";
};

template <typename PreconditionerType, typename OptionsGroup>
struct SubdomainPreconditioner {
  // static std::string name() noexcept {
  //   return option_name<PreconditionerType>();
  // }
  using type = PreconditionerType;
  using group = OptionsGroup;  // SubdomainPreconditionerGroup<OptionsGroup>;
  static constexpr Options::String help =
      "Options for the subdomain solver's preconditioner";
};

}  // namespace OptionTags

/// Tags related to the Schwarz solver
namespace Tags {

/// Number of points a subdomain can overlap with its neighbor
template <typename OptionsGroup>
struct MaxOverlap : db::SimpleTag {
  static std::string name() noexcept {
    return "MaxOverlap(" + Options::name<OptionsGroup>() + ")";
  }
  using type = size_t;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::MaxOverlap<OptionsGroup>>;
  static type create_from_options(const type& value) noexcept { return value; }
};

/// The serial linear solver used to solve subdomain operators
template <typename OptionsGroup>
struct SubdomainSolverBase : db::BaseTag {
  static std::string name() noexcept {
    return "SubdomainSolver(" + Options::name<OptionsGroup>() + ")";
  }
};

/// The serial linear solver of type `SolverType` used to solve subdomain
/// operators
template <typename SolverType, typename OptionsGroup>
struct SubdomainSolver : SubdomainSolverBase<OptionsGroup>, db::SimpleTag {
  using type = SolverType;
  static constexpr bool pass_metavariables = false;
  using option_tags =
      tmpl::list<OptionTags::SubdomainSolver<SolverType, OptionsGroup>>;
  static type create_from_options(const type& value) noexcept { return value; }
};

template <typename OptionsGroup>
struct SubdomainPreconditionerBase : db::BaseTag {};

template <typename PreconditionerType, typename OptionsGroup>
struct SubdomainPreconditioner : db::SimpleTag,
                                 SubdomainPreconditionerBase<OptionsGroup> {
  using type = PreconditionerType;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<
      OptionTags::SubdomainPreconditioner<PreconditionerType, OptionsGroup>>;
  static type create_from_options(const type& value) noexcept { return value; }
};

/*!
 * \brief The `Tag` on the overlap region with each neighbor, i.e. on a region
 * extruding from the central element.
 *
 * Note that data on an overlap with a neighbor is typically oriented according
 * to the neighbor's orientation, so re-orientation needs to happen whenever
 * the data cross element boundaries.
 */
template <typename Tag, size_t Dim, typename OptionsGroup>
struct Overlaps : db::SimpleTag {
  static std::string name() noexcept {
    return "Overlaps(" + db::tag_name<Tag>() + ", " +
           Options::name<OptionsGroup>() + ")";
  }
  using tag = Tag;
  using type = OverlapMap<Dim, typename Tag::type>;
};

/// The number of points a neighbor's subdomain extends into the element
template <size_t Dim, typename OptionsGroup>
struct IntrudingExtents : db::SimpleTag {
  static std::string name() noexcept {
    return "IntrudingExtents(" + Options::name<OptionsGroup>() + ")";
  }
  using type = std::array<size_t, Dim>;
};

/// The width in element-logical coordinates that a neighbor's subdomain extends
/// into the element
template <size_t Dim, typename OptionsGroup>
struct IntrudingOverlapWidths : db::SimpleTag {
  static std::string name() noexcept {
    return "IntrudingOverlapWidths(" + Options::name<OptionsGroup>() + ")";
  }
  using type = std::array<double, Dim>;
};

/// Weighting field for combining data from multiple overlapping subdomains
template <typename OptionsGroup>
struct Weight : db::SimpleTag {
  static std::string name() noexcept {
    return "Weight(" + Options::name<OptionsGroup>() + ")";
  }
  using type = Scalar<DataVector>;
};

/*!
 * \brief A diagnostic quantity to check that weights are conserved
 *
 * This quantity and the `Tags::Weight` on the element should sum to one on all
 * grid points. Residual values indicate that overlap data from neighboring
 * subdomains and data on the element are combined in a non-conservative way.
 */
template <typename OptionsGroup>
struct SummedIntrudingOverlapWeights : db::SimpleTag {
  static std::string name() noexcept {
    return "SummedIntrudingOverlapWeights(" + Options::name<OptionsGroup>() +
           ")";
  }
  using type = Scalar<DataVector>;
};

}  // namespace Tags
}  // namespace LinearSolver::Schwarz

namespace db {
namespace detail {
// This implementation mirrors the interface tag subitems in `Domain/Tags.hpp`.
// Please see that implementation for details.
template <typename VariablesTag, size_t Dim, typename OptionsGroup,
          typename Tags = typename VariablesTag::type::tags_list>
struct OverlapSubitemsImpl;

template <typename VariablesTag, size_t Dim, typename OptionsGroup,
          typename... Tags>
struct OverlapSubitemsImpl<VariablesTag, Dim, OptionsGroup,
                           tmpl::list<Tags...>> {
  using type = tmpl::list<
      LinearSolver::Schwarz::Tags::Overlaps<Tags, Dim, OptionsGroup>...>;
  using tag =
      LinearSolver::Schwarz::Tags::Overlaps<VariablesTag, Dim, OptionsGroup>;
  using return_type = NoSuchType;
  template <typename Subtag>
  static void create_item(
      const gsl::not_null<typename tag::type*> parent_value,
      const gsl::not_null<typename Subtag::type*> sub_value) noexcept {
    sub_value->clear();
    for (auto& [overlap_id, parent_data] : *parent_value) {
      auto& parent_vars = get<typename Subtag::tag>(parent_data);
      auto& sub_var = (*sub_value)[overlap_id];
      for (auto vars_it = parent_vars.begin(), sub_var_it = sub_var.begin();
           vars_it != parent_vars.end(); ++vars_it, ++sub_var_it) {
        sub_var_it->set_data_ref(&*vars_it);
      }
    }
  }
  template <typename Subtag>
  static void create_compute_item(
      const gsl::not_null<typename Subtag::type*> sub_value,
      const typename tag::type& parent_value) noexcept {
    for (const auto& [direction, parent_data] : parent_value) {
      const auto& parent_vars = get<typename Subtag::tag>(parent_data);
      auto& sub_var = (*sub_value)[direction];
      auto sub_var_it = sub_var.begin();
      for (auto vars_it = parent_vars.begin(); vars_it != parent_vars.end();
           ++vars_it, ++sub_var_it) {
        sub_var_it->set_data_ref(const_cast<DataVector*>(&*vars_it));  // NOLINT
      }
    }
  }
};
}  // namespace detail

template <typename VariablesTag, size_t Dim, typename OptionsGroup>
struct Subitems<
    LinearSolver::Schwarz::Tags::Overlaps<VariablesTag, Dim, OptionsGroup>,
    Requires<tt::is_a_v<Variables, typename VariablesTag::type>>>
    : detail::OverlapSubitemsImpl<VariablesTag, Dim, OptionsGroup> {};

}  // namespace db
