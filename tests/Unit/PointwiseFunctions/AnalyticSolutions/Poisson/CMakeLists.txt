# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_PoissonSolutions")

set(LIBRARY_SOURCES
  Test_Lorentzian.cpp
  Test_Moustache.cpp
  Test_ProductOfSinusoids.cpp
  Test_Zero.cpp
  )

add_test_library(
  ${LIBRARY}
  "PointwiseFunctions/AnalyticSolutions/Poisson/"
  "${LIBRARY_SOURCES}"
  "CoordinateMaps;DataStructures;Domain;Poisson\
;PoissonSolutions;Spectral;Utilities"
  )
