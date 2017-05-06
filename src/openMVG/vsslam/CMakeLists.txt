
find_library(SLAMPP_BASE_LIBRARY slampp_base HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CAMD_LIBRARY camd_dlong HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CSPARSE_LIBRARY csparse HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CXSPARSE_LIBRARY cxsparse HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CHOLMOD_LIBRARY cholmod HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_EIGEN_LIBRARY eigen HINTS ${SLAMPP_SOURCE_DIR}/lib)

MESSAGE( STATUS "ceres: ${PROJECT_SOURCE_DIR}/third_party/ceres-solver/include" )

include_directories("${PROJECT_SOURCE_DIR}/third_party/ceres-solver/include")
include_directories("${PROJECT_SOURCE_DIR}/third_party/ceres-solver/config") # !!
include_directories("${PROJECT_SOURCE_DIR}/third_party/ceres-solver/internal/ceres/miniglog") # !!
# swine

IF(NOT OpenMVG_BUILD_OPENGL_EXAMPLES)
	ADD_DEFINITIONS(-DSWINE_NOGL)
ENDIF(NOT OpenMVG_BUILD_OPENGL_EXAMPLES)
# swine

file(
  GLOB_RECURSE
  vsslam_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  vsslam_files_cpp
  *.cpp
)


MESSAGE( STATUS "aa: ${SLAMPP_INCLUDE_DIR}")


MESSAGE( STATUS "qt: ${Qt5Core_DIR}")

add_library(openMVG_vsslam ${vsslam_files_header}  ${SLAMPP_HEADERS} ${vsslam_files_cpp})
target_link_libraries(
  openMVG_vsslam
  openMVG_sfm
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES}
  glog
  openMVG_lInftyComputerVision
  openMVG_system
  openMVG_matching
  vlsift
  ${SLAMPP_BASE_LIBRARY}
  ${SLAMPP_CAMD_LIBRARY}
  ${SLAMPP_CSPARSE_LIBRARY}
  ${SLAMPP_CXSPARSE_LIBRARY}
  ${SLAMPP_CHOLMOD_LIBRARY}
)

set_target_properties(openMVG_vsslam PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_vsslam DESTINATION lib EXPORT openMVG-targets)



# debug - see the properties that vsslam is compiled with
get_directory_property( DirDefs COMPILE_DEFINITIONS )
message( "COMPILE_DEFINITIONS = ${DirDefs}" )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message( "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}" )
	message( "CMAKE_C_FLAGS_DEBUG = ${CMAKE_C_FLAGS_DEBUG}" )
ELSE(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message( "CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}" )
	message( "CMAKE_C_FLAGS_RELEASE = ${CMAKE_C_FLAGS_RELEASE}" )
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
message( "CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}" )
