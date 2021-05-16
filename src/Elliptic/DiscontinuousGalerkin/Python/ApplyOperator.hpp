// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <pybind11/pybind11.h>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Variables.hpp"
#include "Elliptic/DiscontinuousGalerkin/DgOperator.hpp"

namespace py = pybind11;

namespace elliptic::dg::py_bindings {

/// Prefix tag that represents the elliptic DG operator applied to fields.
template <typename Tag>
struct DgOperatorAppliedTo : db::PrefixTag {
  using type = typename Tag::type;
  using tag = Tag;
};

/// An element in a DG domain
template <size_t Dim>
struct DgElement {
  Mesh<Dim> mesh;
  Element<Dim> element;
  ElementMap<Dim, Frame::Inertial> element_map;
  InverseJacobian<DataVector, Dim, Frame::Logical, Frame::Inertial>
      inv_jacobian;
};

/// Defines an ordering of elements by block ID first, then by segment
/// refinement level and last by segment index in each dimension in turn.
template <size_t Dim>
struct ElementOrdering {
  bool operator()(const ElementId<Dim>& lhs,
                  const ElementId<Dim>& rhs) const noexcept;
};

/// An ordered map of `DgElement`s
template <size_t Dim>
using DgElementArray =
    std::map<ElementId<Dim>, DgElement<Dim>, ElementOrdering<Dim>>;

/// Construct a `DgElementArray` from the `domain_creator
template <size_t Dim>
DgElementArray<Dim> create_elements(
    const DomainCreator<Dim>& domain_creator) noexcept;

/// Construct all mortars for the given `element_id`
template <size_t VolumeDim>
::dg::MortarMap<VolumeDim,
                std::pair<Mesh<VolumeDim - 1>, ::dg::MortarSize<VolumeDim - 1>>>
create_mortars(const ElementId<VolumeDim>& element_id,
               const DgElementArray<VolumeDim>& dg_elements) noexcept;

template <typename System>
struct Workspace {
 private:
  static constexpr size_t Dim = System::volume_dim;

 public:
  Variables<typename System::auxiliary_fields> auxiliary_vars;
  Variables<typename System::auxiliary_fluxes> auxiliary_fluxes;
  Variables<typename System::primal_fluxes> primal_fluxes;
  ::dg::MortarMap<
      Dim, ::elliptic::dg::MortarData<size_t, typename System::primal_fields,
                                      typename System::primal_fluxes>>
      all_mortar_data;
};

template <typename System, size_t Dim = System::volume_dim>
std::unordered_map<ElementId<Dim>,
                   Variables<db::wrap_tags_in<DgOperatorAppliedTo,
                                              typename System::primal_fields>>>
apply_operator(
    const gsl::not_null<std::unordered_map<ElementId<Dim>, Workspace<System>>*>
        workspace,
    const DgElementArray<Dim>& dg_element_array,
    const std::unordered_map<
        ElementId<Dim>, Variables<typename System::primal_fields>>& operand,
    const double penalty_parameter) {
  constexpr size_t temporal_id = 0;
  std::unordered_map<ElementId<Dim>,
                     Variables<db::wrap_tags_in<
                         DgOperatorAppliedTo, typename System::primal_fields>>>
      result{};
  for (const auto& [element_id, dg_element] : dg_element_array) {
    auto& this_workspace = workspace->at(element_id);
    const auto& vars = operand.at(element_id);
    ::dg::prepare_mortar_data(
        make_not_null(&this_workspace.auxiliary_vars),
        make_not_null(&this_workspace.auxiliary_fluxes),
        make_not_null(&this_workspace.primal_fluxes),
        make_not_null(&this_workspace.all_mortar_data), vars,
        dg_element.element, dg_element.mesh, dg_element.inv_jacobian,
        dg_element.internal_face_normals, dg_element.external_face_normals,
        dg_element.internal_face_normal_magnitudes,
        dg_element.external_face_normal_magnitudes,
        dg_element.all_mortar_meshes, dg_element.all_mortar_sizes, temporal_id,
        apply_boundary_condition, fluxes_args, sources_args,
        fluxes_args_on_internal_faces, fluxes_args_on_external_faces);
  }
  for (const auto& [element_id, dg_element] : dg_element_array) {
    auto& this_workspace = workspace->at(element_id);
    const auto& vars = operand.at(element_id);
    ::dg::apply_operator(
        make_not_null(&this_result),
        make_not_null(&this_workspace.all_mortar_data), vars,
        this_workspace.primal_fluxes, dg_element.mesh, dg_element.inv_jacobian,
        dg_element.internal_face_normal_magnitudes,
        dg_element.external_face_normal_magnitudes,
        dg_element.all_mortar_meshes, dg_element.all_mortar_sizes,
        penalty_parameter, temporal_id, sources_args);
  }
}

template <typename System>
void bind_apply_operator(py::module& m,
                         const std::string& function_name) {  // NOLINT
  m.def(function_name.c_str(), &apply_operator<System>);
}

}  // namespace elliptic::dg::py_bindings
