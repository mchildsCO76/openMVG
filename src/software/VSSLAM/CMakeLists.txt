
include_directories("${PROJECT_SOURCE_DIR}/third_party/ceres-solver/include")
include_directories("${PROJECT_SOURCE_DIR}/third_party/ceres-solver/config") # !!
# swine

if (OpenMVG_BUILD_OPENGL_EXAMPLES)

  #
  # OpenGL
  #
  find_package(OpenGL REQUIRED)
  include_directories(${OPENGL_INCLUDE_DIR})
  add_definitions(-DGL_GLEXT_PROTOTYPES -DGL3_PROTOTYPES)

  if(APPLE)
    add_definitions(-DGLFW_INCLUDE_GLCOREARB)
  endif()

  #
  # GLFW
  #
  if(!APPLE)
    add_definitions(-DGLFW_NO_GLU)
  endif()

  find_file(GL3HEADER gl3.h)
  if(NOT ${GL3HEADER} MATCHES NOTFOUND)
    add_definitions(-DGLFW_INCLUDE_GL3)
  endif()

  #
  # END - OpenGL
  #

  # - VSSLAM (WIP)
  #


  add_executable(openMVG_main_VSSLAM main_VSSLAM.cpp CGlWindow.hpp ${SLAMPP_HEADERS})
  
  target_link_libraries(openMVG_main_VSSLAM
    ${OPENGL_gl_LIBRARY}
    glfw
    ${GLFW_LIBRARIES}
    ${OpenMVG_LIBRARIES}
  )

  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VSSLAM
        ${OpenCV_LIBS})
  endif(OpenMVG_USE_OPENCV)

  set_property(TARGET openMVG_main_VSSLAM PROPERTY FOLDER OpenMVG/software)



FIND_PACKAGE(Sophus REQUIRED)
include_directories(${Sophus_INCLUDE_DIRS})

add_executable(openMVG_main_VSSLAM_test_sim3 tests/main_VSSLAM_unitTest.cpp ${SLAMPP_HEADERS})
  
  target_link_libraries(openMVG_main_VSSLAM_test_sim3
    ${OPENGL_gl_LIBRARY}
    glfw
    ${GLFW_LIBRARIES}
    ${OpenMVG_LIBRARIES}
    ${Sophus_LIBRARIES}
  )

  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VSSLAM_test_sim3
        ${OpenCV_LIBS})
  endif(OpenMVG_USE_OPENCV)

  set_property(TARGET openMVG_main_VSSLAM_test_sim3 PROPERTY FOLDER OpenMVG/software)

else() # swine

	ADD_DEFINITIONS(SWINE_NOGL)
	# no opengl

  # - VSSLAM (WIP)
  #


  add_executable(openMVG_main_VSSLAM main_VSSLAM.cpp ${SLAMPP_HEADERS}) # CGlWindow.hpp removed
  
  target_link_libraries(openMVG_main_VSSLAM
    ${OpenMVG_LIBRARIES}
  )

	#     ${OPENGL_gl_LIBRARY}     glfw     ${GLFW_LIBRARIES} # removed

  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VSSLAM
        ${OpenCV_LIBS})
  endif(OpenMVG_USE_OPENCV)

  set_property(TARGET openMVG_main_VSSLAM PROPERTY FOLDER OpenMVG/software)



FIND_PACKAGE(Sophus REQUIRED)
include_directories(${Sophus_INCLUDE_DIRS}) # swine - dont have but I put it to path so that it finds it

add_executable(openMVG_main_VSSLAM_test_sim3 tests/main_VSSLAM_unitTest.cpp ${SLAMPP_HEADERS})
  
  target_link_libraries(openMVG_main_VSSLAM_test_sim3
    ${OpenMVG_LIBRARIES}
    ${Sophus_LIBRARIES}
   )

	#     ${OPENGL_gl_LIBRARY}     glfw     ${GLFW_LIBRARIES}   ${Sophus_LIBRARIES} # removed

  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VSSLAM_test_sim3
        ${OpenCV_LIBS})
  endif(OpenMVG_USE_OPENCV)

  set_property(TARGET openMVG_main_VSSLAM_test_sim3 PROPERTY FOLDER OpenMVG/software)


endif()
