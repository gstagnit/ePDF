set(kernel_source
  specialfunctions.cc
  anomalousdimensions.cc
  xdistributions.cc
  ndistributions.cc
  alphaem.cc
  analyticsolutions.cc
  numericintegrals.cc
  )

add_library(kernel OBJECT ${kernel_source})

option(SHARED "Build shared-libray instead of static-libray" ON)
if(SHARED)
  add_library(ePDF SHARED $<TARGET_OBJECTS:kernel>)
else(SHARED)
  add_library(ePDF STATIC $<TARGET_OBJECTS:kernel>)
endif(SHARED)

target_link_libraries(ePDF ${YAML_LDFLAGS} ${GSL_LIBRARIES})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/inc/ePDF DESTINATION include)
install(TARGETS ePDF DESTINATION lib)
