
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

endif()
