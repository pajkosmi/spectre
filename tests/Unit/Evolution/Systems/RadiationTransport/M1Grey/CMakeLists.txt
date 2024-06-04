# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_M1Grey")

set(LIBRARY_SOURCES
  BoundaryConditions/Test_DirichletAnalytic.cpp # success after fixing [s] bug
  BoundaryConditions/Test_Periodic.cpp # success
  BoundaryCorrections/Test_Rusanov.cpp # success
  Test_Actions.cpp # success
  Test_Fluxes.cpp # success
  Test_M1Closure.cpp # 1 failed from broken & untested function
  Test_M1HydroCoupling.cpp # success (will need to be updated for coupling)
  Test_Sources.cpp #pass when including source_n and source_i
  Test_System.cpp # 1 failure with broken untested M1Closure fxn
  Test_Tags.cpp #success
  Test_TimeDerivativeTerms.cpp # passed with fixed sources
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DataStructures
  DataStructuresHelpers
  DiscontinuousGalerkin
  Domain
  M1Grey
  M1GreySolutions
  Parallel
  Spectral
  Utilities
  )
