# get_version : get version information form files
#
# version.txt - version string
# version_date.txt - version date sting
#
# get_version(VER VER_DEV DATE)
#
# will store version string in "VER" and version date string in "DATE"
function(get_version ver_short ver_full ver_date)
  # check version file
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
    file(STRINGS "version.txt" ver_file LIMIT_COUNT 1)
  else()
    set(ver_file "0.0.0-dev0") # default version
  endif()

  # version should be of the form (semver.org):
  # - 1.2.3-dev0 (development with number)
  # - 1.2.3-rc1  (release candidate with number)
  # - 1.2.3      (release)
  # remove possible "v" prefix and find major.minor.patch version
  string(REGEX MATCH "^v?([0-9]+)" _ ${ver_file})
  set(ver_major ${CMAKE_MATCH_1})
  string(REGEX MATCH "^v?[0-9]+\.([0-9]+)" _ ${ver_file})
  set(ver_minor ${CMAKE_MATCH_1})
  string(REGEX MATCH "^v?[0-9]+\.[0-9]+\.([0-9]+)" _ ${ver_file})
  set(ver_patch ${CMAKE_MATCH_1})
  # find pre-release tag
  string(REGEX MATCH ".*-(.+)" _ ${ver_file})
  set(ver_pre ${CMAKE_MATCH_1})

  # create the version string for cmake (fill up with 0)
  if ("${ver_major}" STREQUAL "")
    set(ver_major 0) # default version
  endif()
  if ("${ver_minor}" STREQUAL "")
    set(ver_minor 0) # default version
  endif()
  if ("${ver_patch}" STREQUAL "")
    set(ver_patch 0) # default version
  endif()
  # create version x.y.z
  set(ver ${ver_major}.${ver_minor}.${ver_patch})

  # check date file (if not a development version)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version_date.txt" AND (NOT (${ver_pre} MATCHES "^dev.*")))
    file(STRINGS "version_date.txt" date LIMIT_COUNT 1)
  else()
    string(TIMESTAMP date "%Y-%m-%d") # current date
  endif()

  message(STATUS "use VERSION: ${ver} (from ${ver_file})")
  message(STATUS "use DATE:    ${date}")

  # set given variables in parent scope
  set(${ver_date} ${date} PARENT_SCOPE)
  set(${ver_full} ${ver_file} PARENT_SCOPE)
  set(${ver_short} ${ver} PARENT_SCOPE)
endfunction()
