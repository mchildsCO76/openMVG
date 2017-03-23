
IF (OpenMVG_BUILD_OPENGL_EXAMPLES)

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

set(SLAMPP_SOURCE_DIR "/home/klemen/workspace/ANU/SlamPP_Dev/slam-plus-plus-dev_code/trunk/SLAM_plus_plus")
set(SLAMPP_INCLUDE_DIR ${SLAMPP_SOURCE_DIR}/include ${SLAMPP_SOURCE_DIR}/include/cholmod/SuiteSparse)

include_directories(${SLAMPP_INCLUDE_DIR})


#file(GLOB_RECURSE SLAMPP_HEADERS ${SLAMPP_SOURCE_DIR}/include/slam/*.h)

  add_executable(openMVG_main_VSSLAM main_VSSLAM.cpp ../VO/CGlWindow.hpp ${SLAMPP_HEADERS})
  target_link_libraries(openMVG_main_VSSLAM
    ${OPENGL_gl_LIBRARY}
    glfw
    ${GLFW_LIBRARIES}
    ${OpenMVG_LIBRARIES}
    vlsift
  )

  IF(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VSSLAM
        ${OpenCV_LIBS})
  ENDIF(OpenMVG_USE_OPENCV)

  SET_PROPERTY(TARGET openMVG_main_VSSLAM PROPERTY FOLDER OpenMVG/software)

ENDIF()
