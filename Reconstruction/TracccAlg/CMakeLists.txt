# Headers and Libraries

# Modules
gaudi_add_module(TracccAlg
                 SOURCES src/VXDPredict/VXDPredict.cpp
                 LINK DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiAlgLib Gaudi::GaudiKernel
                      ${LCIO_LIBRARIES} 
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${podio_LIBRARIES} podio::podioRootIO
                      /home/zhangyz/CEPCSW/Reconstruction/TracccAlg/lib/libsyclevent.so
)

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
