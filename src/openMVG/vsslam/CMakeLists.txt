
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


add_library(openMVG_vsslam ${vsslam_files_header} ${vsslam_files_cpp})
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
)
set_target_properties(openMVG_vsslam PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_vsslam DESTINATION lib EXPORT openMVG-targets)


add_subdirectory(detection)
add_subdirectory(mapping)
add_subdirectory(optimization)
add_subdirectory(tracking)
add_subdirectory(matching)
