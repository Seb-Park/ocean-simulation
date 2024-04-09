# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles/StaticGLEW_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/StaticGLEW_autogen.dir/ParseCache.txt"
  "CMakeFiles/ocean_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/ocean_autogen.dir/ParseCache.txt"
  "StaticGLEW_autogen"
  "ocean_autogen"
  )
endif()
