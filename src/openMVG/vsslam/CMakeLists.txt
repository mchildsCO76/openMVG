
set(SLAMPP_SOURCE_DIR "/home/klemen/workspace/ANU/SlamPP_Dev/slam-plus-plus-dev_code/trunk/SLAM_plus_plus")
set(SLAMPP_INCLUDE_DIR ${SLAMPP_SOURCE_DIR}/include ${SLAMPP_SOURCE_DIR}/include/cholmod/SuiteSparse)

include_directories(${SLAMPP_INCLUDE_DIR})
find_library(SLAMPP_BASE_LIBRARY slampp_base HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CAMD_LIBRARY camd_dlong HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CSPARSE_LIBRARY csparse HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CXSPARSE_LIBRARY cxsparse HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_CHOLMOD_LIBRARY cholmod HINTS ${SLAMPP_SOURCE_DIR}/lib)
find_library(SLAMPP_EIGEN_LIBRARY eigen HINTS ${SLAMPP_SOURCE_DIR}/lib)


file(GLOB_RECURSE SLAMPP_HEADERS ${SLAMPP_SOURCE_DIR}/include/slam/*.h)
set(SLAM_CAMD_LIBRARY /home/klemen/workspace/ANU/SlamPP_Dev/slam-plus-plus-dev_code/trunk/SLAM_plus_plus/lib/libcamd_dlong.a)

add_subdirectory(detection)
add_subdirectory(mapping)
add_subdirectory(optimization)
add_subdirectory(tracking)
add_subdirectory(matching)

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

MESSAGE( STATUS "aa: ${SLAMPP_CAMD_LIBRARY}")
add_library(openMVG_vsslam ${vsslam_files_header} ${SLAMPP_HEADERS} ${vsslam_files_cpp}  )

target_link_libraries(
  openMVG_vsslam
  openMVG_sfm
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES}
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


