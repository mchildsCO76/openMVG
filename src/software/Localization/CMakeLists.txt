
###
# Localize images in an existing SfM reconstruction
###
ADD_EXECUTABLE(openMVG_main_SfM_Localization main_SfM_Localization.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_Localization
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )
  

ADD_EXECUTABLE(openMVG_main_SfM_ReconstructionLocalization main_SfM_ReconstructionLocalization.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_ReconstructionLocalization
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )
  
# Installation rules
SET_PROPERTY(TARGET openMVG_main_SfM_Localization PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_Localization DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_SfM_ReconstructionLocalization PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_ReconstructionLocalization DESTINATION bin/)
