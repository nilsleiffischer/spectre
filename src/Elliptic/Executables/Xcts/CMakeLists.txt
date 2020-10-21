# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  Elasticity
  Informer
  IO
  LinearOperators
  Options
  ParallelLinearSolver
  ParallelMultigrid
  ParallelNonlinearSolver
  ParallelSchwarz
  Poisson
  Utilities
  Xcts
  XctsBoundaryConditions
  )

function(add_xcts_executable
  EXECUTABLE_NAME XCTS_SYSTEM BACKGROUND BOUNDARY_CONDITIONS INITIAL_GUESS)
  add_spectre_parallel_executable(
    "Solve${EXECUTABLE_NAME}"
    SolveXctsProblem
    Elliptic/Executables/Xcts
    "Metavariables<${XCTS_SYSTEM}, ${BACKGROUND}, ${BOUNDARY_CONDITIONS}, ${INITIAL_GUESS}>"
    "${LIBS_TO_LINK}"
    )
endfunction(add_xcts_executable)

# Poisson correction system with ProductOfSinusoids solution
add_xcts_executable(
  PoissonCorrection1D
  "Poisson::FirstOrderSystem<1, Poisson::Geometry::Euclidean>"
  Poisson::Solutions::ProductOfSinusoids<1>
  Poisson::Solutions::ProductOfSinusoids<1>
  Poisson::Solutions::Zero<1>
  )
target_link_libraries(
  SolvePoissonCorrection1D
  PRIVATE
  PoissonSolutions
  )

# Xcts system with ConstantDensityStar solution
add_xcts_executable(
  ConstantDensityStar
  "Xcts::FirstOrderSystem<Xcts::Equations::Hamiltonian, Xcts::Geometry::Euclidean>"
  Xcts::Solutions::ConstantDensityStar
  Xcts::Solutions::ConstantDensityStar
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveConstantDensityStar
  PRIVATE
  XctsSolutions
  )

# Xcts system with single black-hole solutions in various coordinates
# (i) Conformally-flat maximally-sliced Schwarzschild (works for
# Hamiltonian&Lapse equations, broken when including momentum constraint and
# AH boundary conditions because shift condition conflicts with beta=0 solution)
add_xcts_executable(
  SchwarzschildIsotropic
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::Euclidean>"
  Xcts::Solutions::Schwarzschild<Xcts::Solutions::SchwarzschildCoordinates::Isotropic>
  Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::Euclidean>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveSchwarzschildIsotropic
  PRIVATE
  XctsSolutions
  )
# (ii) Conformally-flat Schwarzschild with Kerr-Schild slicing
add_xcts_executable(
  SchwarzschildKSI
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::Euclidean>"
  Xcts::Solutions::Schwarzschild<Xcts::Solutions::SchwarzschildCoordinates::KerrSchildIsotropic>
  Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::Euclidean>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveSchwarzschildKSI
  PRIVATE
  XctsSolutions
  )
# (iii) Kerr solution in Kerr-Schild coordinates (BROKEN: perhaps non-Euclidean
# term in AH b.c. is wrong)
add_xcts_executable(
  KerrSchild
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::NonEuclidean>"
  Xcts::Solutions::Kerr<Xcts::Solutions::KerrCoordinates::KerrSchild>
  # Xcts::Solutions::Kerr<Xcts::Solutions::KerrCoordinates::KerrSchild>
  Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::NonEuclidean>
  # Xcts::Solutions::Flatness
  Xcts::Solutions::Kerr<Xcts::Solutions::KerrCoordinates::KerrSchild>
  )
target_link_libraries(
  SolveKerrSchild
  PRIVATE
  XctsSolutions
  )

# BBH initial data with various choices for the free data
# (i) Conformally flat and maximally sliced (BROKEN: see single BH (i))
add_xcts_executable(
  BbhCFMS
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::Euclidean>"
  Xcts::AnalyticData::BlackHoleBinary<Xcts::AnalyticData::BackgroundSpacetime::FlatMaximallySliced>
  Xcts::BoundaryConditions::Binary<Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::Euclidean>>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveBbhCFMS
  PRIVATE
  XctsAnalyticData
  )
# (ii) Conformally flat with Eddington-Finkelstein slicing
add_xcts_executable(
  BbhKSI
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::Euclidean>"
  Xcts::AnalyticData::BlackHoleBinary<Xcts::AnalyticData::BackgroundSpacetime::FlatKerrSchildIsotropicSliced>
  Xcts::BoundaryConditions::Binary<Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::Euclidean>>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveBbhKSI
  PRIVATE
  XctsAnalyticData
  )
# (iii) Superposed Kerr-Schild (TODO: Handle non-Euclidean magnitude in
# subdomain operator)
add_xcts_executable(
  BbhSKS
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::NonEuclidean>"
  Xcts::AnalyticData::BlackHoleBinary<Xcts::AnalyticData::BackgroundSpacetime::SuperposedKerrSchild>
  Xcts::BoundaryConditions::Binary<Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::NonEuclidean>>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveBbhSKS
  PRIVATE
  XctsAnalyticData
  )
# (iv) Superposed Harmonic Kerr (TODO: Implement Kerr solution in harmonic
# coords)
add_xcts_executable(
  BbhSHK
  "Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift, Xcts::Geometry::NonEuclidean>"
  Xcts::AnalyticData::BlackHoleBinary<Xcts::AnalyticData::BackgroundSpacetime::SuperposedHarmonic>
  Xcts::BoundaryConditions::Binary<Xcts::BoundaryConditions::ApparentHorizon<Xcts::Geometry::NonEuclidean>>
  Xcts::Solutions::Flatness
  )
target_link_libraries(
  SolveBbhSHK
  PRIVATE
  XctsAnalyticData
  )