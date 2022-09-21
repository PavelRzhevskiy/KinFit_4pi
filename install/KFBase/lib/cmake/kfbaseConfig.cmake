
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was kfbaseConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

####################################################################################

find_package(Eigen3 REQUIRED NO_MODULE)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
include("${CMAKE_CURRENT_LIST_DIR}/kfbaseTargets.cmake" )
set(KFBASE_ROOT_DIR /spoolA/przhevskiy/kinfit/install/KFBase)
set(KFBASE_INCLUDE_DIRS ${KFBASE_ROOT_DIR}/include)
set(KFBASE_LIBRARIES
  ${KFBASE_ROOT_DIR}/lib/libkfbase_newtonian_opt.so
  ${KFBASE_ROOT_DIR}/lib/libkfbase_core.so)
