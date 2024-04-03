# Headers and Libraries

# Set up externs
# find_package( vecmem 0.7.0 REQUIRED COMPONENTS LANGUAGE )
# find_package( Eigen3 REQUIRED )
# find_package( Thrust REQUIRED )

# # Set up an IMPORTED library on top of the Thrust library/libraries. One that
# # the TRACCC/Detray code could depend on publicly.
# set( TRACCC_THRUST_OPTIONS "" CACHE STRING
#    "Extra options for configuring how Thrust should be used" )
# mark_as_advanced( TRACCC_THRUST_OPTIONS )
# thrust_create_target( traccc::Thrust ${TRACCC_THRUST_OPTIONS} )
# # Make sure that Thrust/CUB headers are treated as system headers.
# get_target_property( _thrustIncludeDir Thrust::Thrust INTERFACE_INCLUDE_DIRECTORIES )
# get_target_property( _cubIncludeDir CUB::CUB INTERFACE_INCLUDE_DIRECTORIES )
# target_include_directories( traccc::Thrust
#    SYSTEM INTERFACE ${_thrustIncludeDir} ${_cubIncludeDir} )
# unset( _thrustIncludeDir )
# unset( _cubIncludeDir )

# find_package( TBB REQUIRED )
# find_package( Kokkos REQUIRED )
# find_package( alpaka REQUIRED )
# find_package( algebra-plugins REQUIRED )
# find_package( covfie REQUIRED )
# find_package( dfelibs REQUIRED )
# find_package( detray REQUIRED )
# find_package( Acts REQUIRED )
# find_package( GTest REQUIRED )

enable_language( SYCL )
# Modules
gaudi_add_module(TracccAlg
                 SOURCES src/TracccAlg.cpp
                 LINK DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiAlgLib Gaudi::GaudiKernel
                      ${LCIO_LIBRARIES} 
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${podio_LIBRARIES} podio::podioRootIO
)

target_link_libraries(TracccAlg PUBLIC traccc::options vecmem::core vecmem::sycl traccc::io
                      traccc::core traccc::device_common traccc::sycl traccc::performance)

if (GenFit_FOUND)
  target_link_libraries(TracccAlg PUBLIC GenFit::genfit2)
endif()

target_include_directories(TracccAlg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS TracccAlg
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
