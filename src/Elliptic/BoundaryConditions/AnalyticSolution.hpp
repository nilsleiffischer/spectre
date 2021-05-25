// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <ostream>
#include <pup.h>
#include <string>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FaceNormal.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryConditionType.hpp"
#include "Elliptic/BoundaryConditions/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/NormalDotFlux.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace elliptic::BoundaryConditions {

/// \cond
template <typename System, size_t Dim, typename FieldTags, typename FluxTags,
          typename Registrars>
struct AnalyticSolution;

namespace Registrars {
template <typename System, size_t Dim = System::volume_dim,
          typename FieldTags = typename System::primal_fields,
          typename FluxTags = typename System::primal_fluxes>
struct AnalyticSolution {
  template <typename Registrars>
  using f = BoundaryConditions::AnalyticSolution<System, Dim, FieldTags,
                                                 FluxTags, Registrars>;
};
}  // namespace Registrars

template <typename System, size_t Dim = System::volume_dim,
          typename FieldTags = typename System::primal_fields,
          typename FluxTags = typename System::primal_fluxes,
          typename Registrars =
              tmpl::list<BoundaryConditions::Registrars::AnalyticSolution<
                  System, Dim, FieldTags, FluxTags>>>
struct AnalyticSolution;
/// \endcond

/*!
 * \brief Impose the analytic solution on the boundary. Works only if an
 * analytic solution exists.
 *
 * The analytic solution is retrieved from `::Tags::AnalyticSolutionsBase`. It
 * must hold solutions for both the `System::primal_fields` and the
 * `System::primal_fluxes`. The user can select to impose the analytic solution
 * as Dirichlet or Neumann boundary conditions for each field separately.
 * Dirichlet boundary conditions are imposed on the fields and Neumann boundary
 * conditions are imposed on the fluxes.
 *
 * See `elliptic::Actions::InitializeAnalyticSolutions` for an action that can
 * add the analytic solutions to the DataBox.
 */
template <typename System, size_t Dim, typename... FieldTags,
          typename... FluxTags, typename Registrars>
class AnalyticSolution<System, Dim, tmpl::list<FieldTags...>,
                       tmpl::list<FluxTags...>, Registrars>
    : public BoundaryCondition<Dim, Registrars> {
 private:
  using Base = BoundaryCondition<Dim, Registrars>;

 public:
  using options =
      tmpl::list<elliptic::OptionTags::BoundaryConditionType<FieldTags>...>;
  static constexpr Options::String help =
      "Boundary conditions from the analytic solution";

  AnalyticSolution() = default;
  AnalyticSolution(const AnalyticSolution&) noexcept = default;
  AnalyticSolution& operator=(const AnalyticSolution&) noexcept = default;
  AnalyticSolution(AnalyticSolution&&) noexcept = default;
  AnalyticSolution& operator=(AnalyticSolution&&) noexcept = default;
  ~AnalyticSolution() noexcept = default;

  /// \cond
  explicit AnalyticSolution(CkMigrateMessage* m) noexcept : Base(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(AnalyticSolution);
  /// \endcond

  /// Select which `elliptic::BoundaryConditionType` to apply for each field
  explicit AnalyticSolution(
      // This pack expansion repeats the type `elliptic::BoundaryConditionType`
      // for each system field
      const typename elliptic::OptionTags::BoundaryConditionType<
          FieldTags>::type... boundary_condition_types) noexcept
      : boundary_condition_types_{boundary_condition_types...} {}

  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const noexcept override {
    return std::make_unique<AnalyticSolution>(*this);
  }

  const auto& boundary_condition_types() const noexcept {
    return boundary_condition_types_;
  }

  using argument_tags = tmpl::list<
      ::Tags::AnalyticSolutionsOnBoundaryBase, domain::Tags::Direction<Dim>,
      ::Tags::Normalized<
          domain::Tags::UnnormalizedFaceNormal<Dim, Frame::Inertial>>>;
  using volume_tags = tmpl::list<::Tags::AnalyticSolutionsOnBoundaryBase>;

  template <typename OptionalAnalyticSolutions>
  void apply(const gsl::not_null<typename FieldTags::type*>... fields,
             const gsl::not_null<typename FieldTags::type*>... n_dot_fluxes,
             const OptionalAnalyticSolutions& optional_analytic_solutions,
             const Direction<Dim>& direction,
             const tnsr::i<DataVector, Dim>& face_normal) const noexcept {
    const auto& analytic_solutions = [&optional_analytic_solutions]() noexcept
        -> const auto& {
      if constexpr (tt::is_a_v<std::optional, OptionalAnalyticSolutions>) {
        if (not optional_analytic_solutions.has_value()) {
          ERROR_NO_TRACE(
              "Trying to impose boundary conditions from an analytic solution, "
              "but no analytic solution is available. You probably selected "
              "the 'AnalyticSolution' boundary condition but chose to solve a "
              "problem that has no analytic solution. If this is the case, you "
              "should probably select a different boundary condition.");
        }
        return *optional_analytic_solutions;
      } else {
        return optional_analytic_solutions;
      }
    }
    ();
    const auto impose_boundary_condition = [this, &analytic_solutions,
                                            &direction, &face_normal](
                                               auto field_tag_v,
                                               auto flux_tag_v,
                                               const auto field,
                                               const auto n_dot_flux) noexcept {
      using field_tag = decltype(field_tag_v);
      using flux_tag = decltype(flux_tag_v);
      switch (get<elliptic::Tags::BoundaryConditionType<field_tag>>(
          boundary_condition_types())) {
        case elliptic::BoundaryConditionType::Dirichlet:
          *field = get<::Tags::Analytic<field_tag>>(
              analytic_solutions.at(direction));
          break;
        case elliptic::BoundaryConditionType::Neumann:
          normal_dot_flux(n_dot_flux, face_normal,
                          get<::Tags::Analytic<flux_tag>>(
                              analytic_solutions.at(direction)));
          break;
        default:
          ERROR("Unsupported boundary condition type: "
                << get<elliptic::Tags::BoundaryConditionType<field_tag>>(
                       boundary_condition_types()));
      }
    };
    EXPAND_PACK_LEFT_TO_RIGHT(impose_boundary_condition(FieldTags{}, FluxTags{},
                                                        fields, n_dot_fluxes));
  }

  using argument_tags_linearized = tmpl::list<>;
  using volume_tags_linearized = tmpl::list<>;

  void apply_linearized(
      const gsl::not_null<typename FieldTags::type*>... fields,
      const gsl::not_null<typename FieldTags::type*>... n_dot_fluxes)
      const noexcept {
    const auto impose_boundary_condition = [this](
                                               auto field_tag_v,
                                               const auto field,
                                               const auto n_dot_flux) noexcept {
      using field_tag = decltype(field_tag_v);
      switch (get<elliptic::Tags::BoundaryConditionType<field_tag>>(
          boundary_condition_types())) {
        case elliptic::BoundaryConditionType::Dirichlet:
          for (auto& field_component : *field) {
            field_component = 0.;
          }
          break;
        case elliptic::BoundaryConditionType::Neumann:
          for (auto& n_dot_flux_component : *n_dot_flux) {
            n_dot_flux_component = 0.;
          }
          break;
        default:
          ERROR("Unsupported boundary condition type: "
                << get<elliptic::Tags::BoundaryConditionType<field_tag>>(
                       boundary_condition_types()));
      }
    };
    EXPAND_PACK_LEFT_TO_RIGHT(
        impose_boundary_condition(FieldTags{}, fields, n_dot_fluxes));
  }

  // NOLINTNEXTLINE
  void pup(PUP::er& p) noexcept override;

 private:
  friend bool operator==(const AnalyticSolution& lhs,
                         const AnalyticSolution& rhs) noexcept {
    return lhs.boundary_condition_types_ == rhs.boundary_condition_types_;
  }

  friend bool operator!=(const AnalyticSolution& lhs,
                         const AnalyticSolution& rhs) noexcept {
    return not(lhs == rhs);
  }

  tuples::TaggedTuple<elliptic::Tags::BoundaryConditionType<FieldTags>...>
      boundary_condition_types_{};
};

template <typename System, size_t Dim, typename... FieldTags,
          typename... FluxTags, typename Registrars>
void AnalyticSolution<System, Dim, tmpl::list<FieldTags...>,
                      tmpl::list<FluxTags...>,
                      Registrars>::pup(PUP::er& p) noexcept {
  Base::pup(p);
  p | boundary_condition_types_;
}

/// \cond
template <typename System, size_t Dim, typename... FieldTags,
          typename... FluxTags, typename Registrars>
PUP::able::PUP_ID AnalyticSolution<System, Dim, tmpl::list<FieldTags...>,
                                   tmpl::list<FluxTags...>,
                                   Registrars>::my_PUP_ID = 0;  // NOLINT
/// \endcond

}  // namespace elliptic::BoundaryConditions
