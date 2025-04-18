# Try to find SuperLU 7.x
# Once done, this will define
#
#  SuperLU_FOUND        - system has SuperLU
#  SuperLU_LIBRARY      - Link this to use SuperLU
#  SuperLU_INCLUDE_DIR  - directory of SuperLU headers

find_path(SuperLU_INCLUDE_DIR slu_ddefs.h
  PATHS ${CMAKE_SYSTEM_INCLUDE_PATH} /usr/include /usr/local/include /opt/local/include
  PATH_SUFFIXES superlu SuperLU ""
)

find_library(SuperLU_LIBRARY
  NAMES superlu
  PATHS ${CMAKE_SYSTEM_LIBRARY_PATH} /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib
)

set(SuperLU_FOUND NO)

if(NOT SuperLU_FIND_QUIETLY)
  message(STATUS "Looking for SuperLU version 7")
endif()

if(SuperLU_LIBRARY AND SuperLU_INCLUDE_DIR)
  
  # Now determine if the version is version 7.
  # This voodoo is required as SuperLU developers like to keep things messy.
  
  if(EXISTS "${SuperLU_INCLUDE_DIR}/slu_ddefs.h")
    
    file(READ "${SuperLU_INCLUDE_DIR}/slu_ddefs.h" _slu_HEADER_CONTENTS)
    
    string(REGEX REPLACE ".*version ([0-9]).*" "\\1" SLU_VERSION_MAJOR "${_slu_HEADER_CONTENTS}")
    
    if("${SLU_VERSION_MAJOR}" EQUAL "7")
      set(SuperLU_FOUND YES)
    else()
      if(NOT SuperLU_FIND_QUIETLY)
        message(STATUS "Found an older version of SuperLU, but version 7 is required")
      endif()
    endif()
  else()
    if(NOT SuperLU_FIND_QUIETLY)
      message(STATUS "Found SuperLU and headers, but could not verify version 7")
    endif()
  endif()
endif()

if(SuperLU_LIBRARY AND NOT SuperLU_INCLUDE_DIR)
  if(NOT SuperLU_FIND_QUIETLY)
    message(STATUS "Incomplete installation of SuperLU: no SuperLU headers found")
  endif()
endif()
  
if(NOT SuperLU_FOUND) 
  if(NOT SuperLU_FIND_QUIETLY)
    message(STATUS "Could not find SuperLU")
  endif()
else()
  message(STATUS "Found SuperLU: ${SuperLU_LIBRARY}")
endif()

