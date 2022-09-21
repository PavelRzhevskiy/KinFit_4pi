#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "kfbase::kfbase_newtonian_opt" for configuration "Release"
set_property(TARGET kfbase::kfbase_newtonian_opt APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(kfbase::kfbase_newtonian_opt PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkfbase_newtonian_opt.so"
  IMPORTED_SONAME_RELEASE "libkfbase_newtonian_opt.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS kfbase::kfbase_newtonian_opt )
list(APPEND _IMPORT_CHECK_FILES_FOR_kfbase::kfbase_newtonian_opt "${_IMPORT_PREFIX}/lib/libkfbase_newtonian_opt.so" )

# Import target "kfbase::kfbase_core" for configuration "Release"
set_property(TARGET kfbase::kfbase_core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(kfbase::kfbase_core PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkfbase_core.so"
  IMPORTED_SONAME_RELEASE "libkfbase_core.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS kfbase::kfbase_core )
list(APPEND _IMPORT_CHECK_FILES_FOR_kfbase::kfbase_core "${_IMPORT_PREFIX}/lib/libkfbase_core.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
