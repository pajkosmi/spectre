# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  PRIVATE
  Actions
  Charmxx::main
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  Events
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  GeneralRelativitySolutions
  Hydro
  IO
  Imex
  Informer
  Limiters
  LinearOperators
  MathFunctions
  M1Grey
  M1GreyAnalyticData
  M1GreySolutions
  Observer
  Options
  Parallel
  PhaseControl
  Serialization
  Time
  Utilities
  )

set(EXECUTABLE EvolveM1Grey)

add_spectre_executable(
  ${EXECUTABLE}
  EXCLUDE_FROM_ALL
  EvolveM1Grey.cpp
  )

target_link_libraries(${EXECUTABLE} PRIVATE ${LIBS_TO_LINK})

# add_spectre_executable(
#   EvolveM1GreyConstantM1
#   EvolveM1Grey
#   Evolution/Executables/RadiationTransport/M1Grey
#   EvolutionMetavars
#   "${LIBS_TO_LINK}"
#   )

# add_spectre_executable(
#   EvolveM1GreyHomogeneousSphere
#   EvolveM1Grey
#   Evolution/Executables/RadiationTransport/M1Grey
#   EvolutionMetavars
#   "${LIBS_TO_LINK}"
#   )
