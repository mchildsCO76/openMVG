
###
# Localize images in an existing SfM reconstruction
###
ADD_EXECUTABLE(openMVG_main_SfM_LocalizationImage main_SfM_LocalizationImage.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_LocalizationImage
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )
  

ADD_EXECUTABLE(openMVG_main_SfM_LocalizationReconstruction main_SfM_LocalizationReconstruction.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_LocalizationReconstruction
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )
  
# Installation rules
SET_PROPERTY(TARGET openMVG_main_SfM_LocalizationImage PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_LocalizationImage DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_SfM_LocalizationReconstruction PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_LocalizationReconstruction DESTINATION bin/)
