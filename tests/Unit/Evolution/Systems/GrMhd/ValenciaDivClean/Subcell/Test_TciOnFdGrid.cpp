// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Subcell/TciOnFdGrid.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Utilities/Gsl.hpp"

namespace {
enum class TestThis { AllGood, NeededFixing, PerssonTildeD, PerssonTildeTau };

void test(const TestThis test_this) {
  const Mesh<3> mesh{6, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const double persson_exponent = 5.0;

  auto box = db::create<
      db::AddSimpleTags<evolution::dg::subcell::Tags::Inactive<
                            grmhd::ValenciaDivClean::Tags::TildeD>,
                        evolution::dg::subcell::Tags::Inactive<
                            grmhd::ValenciaDivClean::Tags::TildeTau>,
                        grmhd::ValenciaDivClean::Tags::VariablesNeededFixing,
                        domain::Tags::Mesh<3>>>(
      Scalar<DataVector>(mesh.number_of_grid_points(), 1.0),
      Scalar<DataVector>(mesh.number_of_grid_points(), 1.0),
      test_this == TestThis::NeededFixing, mesh);

  const size_t point_to_change = mesh.number_of_grid_points() / 2;
  if (test_this == TestThis::PerssonTildeTau) {
    db::mutate<evolution::dg::subcell::Tags::Inactive<
        grmhd::ValenciaDivClean::Tags::TildeTau>>(
        make_not_null(&box), [point_to_change](const auto tilde_tau_ptr) {
          get(*tilde_tau_ptr)[point_to_change] *= 2.0;
        });
  } else if (test_this == TestThis::PerssonTildeD) {
    db::mutate<evolution::dg::subcell::Tags::Inactive<
        grmhd::ValenciaDivClean::Tags::TildeD>>(
        make_not_null(&box), [point_to_change](const auto tilde_d_ptr) {
          get(*tilde_d_ptr)[point_to_change] *= 2.0;
        });
  }

  const bool result =
      db::mutate_apply<grmhd::ValenciaDivClean::subcell::TciOnFdGrid>(
          make_not_null(&box), persson_exponent);
  if (test_this == TestThis::AllGood) {
    CHECK_FALSE(result);
  } else {
    CHECK(result);
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.ValenciaDivClean.Subcell.TciOnFdGrid",
                  "[Unit][Evolution]") {
  for (const TestThis& test_this :
       {TestThis::AllGood, TestThis::NeededFixing, TestThis::PerssonTildeD,
        TestThis::PerssonTildeTau}) {
    test(test_this);
  }
}
