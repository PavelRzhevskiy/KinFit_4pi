
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was kfcmdConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

####################################################################################

find_package(kfbase REQUIRED)
include("${CMAKE_CURRENT_LIST_DIR}/kfcmdTargets.cmake" )
set(KFCMD_ROOT_DIR /spoolA/przhevskiy/kinfit/install/KFCmd)
set(KFCMD_INCLUDE_DIRS ${KFCMD_ROOT_DIR}/include)
set(KFCMD_LIBRARIES
  ${KFCMD_ROOT_DIR}/lib/libkfcmd_core.so
  ${KFCMD_ROOT_DIR}/lib/libkfcmd_hypos.so)
