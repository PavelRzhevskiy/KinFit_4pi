#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "kfcmd::kfcmd_core" for configuration "Release"
set_property(TARGET kfcmd::kfcmd_core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(kfcmd::kfcmd_core PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkfcmd_core.so"
  IMPORTED_SONAME_RELEASE "libkfcmd_core.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS kfcmd::kfcmd_core )
list(APPEND _IMPORT_CHECK_FILES_FOR_kfcmd::kfcmd_core "${_IMPORT_PREFIX}/lib/libkfcmd_core.so" )

# Import target "kfcmd::kfcmd_hypos" for configuration "Release"
set_property(TARGET kfcmd::kfcmd_hypos APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(kfcmd::kfcmd_hypos PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkfcmd_hypos.so"
  IMPORTED_SONAME_RELEASE "libkfcmd_hypos.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS kfcmd::kfcmd_hypos )
list(APPEND _IMPORT_CHECK_FILES_FOR_kfcmd::kfcmd_hypos "${_IMPORT_PREFIX}/lib/libkfcmd_hypos.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
