# Visual Studio detection
if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
  OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
  
  set (VS32 TRUE)
  set (VS64 FALSE)
  message (STATUS "Generator: Visual Studio 32 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
  OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
  
if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")
  
  set (VS32 FALSE)
  set (VS64 TRUE)
  message (STATUS "Generator: Visual Studio 64 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")

# Convert environment variables
if (NOT $ENV{LIBRARIES_DIR} STREQUAL "")
	STRING(REGEX REPLACE "\\\\" "/" LIBRARIES_DIR $ENV{LIBRARIES_DIR}) 
endif (NOT $ENV{LIBRARIES_DIR} STREQUAL "")

MACRO(COPY_FILE_IF_CHANGED in_file out_file target)
    IF(${in_file} IS_NEWER_THAN ${out_file})    
        #message("Copying file: ${in_file} to: ${out_file}")
        ADD_CUSTOM_COMMAND (
                TARGET     ${target}
                POST_BUILD
                COMMAND    ${CMAKE_COMMAND}
                ARGS       -E copy ${in_file} ${out_file}
        )
        ENDIF(${in_file} IS_NEWER_THAN ${out_file})
ENDMACRO(COPY_FILE_IF_CHANGED)

MACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED in_file out_dir target)
        GET_FILENAME_COMPONENT(file_name ${in_file} NAME)
        COPY_FILE_IF_CHANGED(${in_file} ${out_dir}/${file_name}
${target})      
ENDMACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED)

#Copies all the files from in_file_list into the out_dir. 
# sub-trees are ignored (files are stored in same out_dir)
MACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED in_file_list out_dir target)
    FOREACH(in_file ${in_file_list})
                COPY_FILE_INTO_DIRECTORY_IF_CHANGED(${in_file}
${out_dir} ${target})
        ENDFOREACH(in_file)     
ENDMACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED)

MACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY in_file target)
	IF (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}/Debug
			${target}
		)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			"${CMAKE_CURRENT_SOURCE_DIR}/${in_file}"
			${EXECUTABLE_OUTPUT_PATH}/Release
			${target}
		)
	ELSE (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}
			${target}
		)
	ENDIF (WIN32)
ENDMACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY)
