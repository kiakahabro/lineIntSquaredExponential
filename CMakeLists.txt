cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Set the project name
# file(STRINGS project.txt PRJ_NAME LIMIT_COUNT 1)
set(PRJ_NAME "double_line_integral")
project(${PRJ_NAME}) 

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Force coloured compiler output
add_compile_options(-fdiagnostics-color)

# CXX flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_CXX_FLAGS_DEBUG " -g -O0")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Compiling for release, disabling debugging")
    add_definitions(-DNDEBUG)
    add_definitions(-DEIGEN_NO_DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Compiling for debugging")
    add_definitions(-DDEBUG)
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
else()
    message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} not recognised.")
endif()



# Vectorisation
# https://eigen.tuxfamily.org/index.php?title=FAQ#How_can_I_enable_vectorization.3F
# Fast math : https://scicomp.stackexchange.com/questions/31216/lapack-equivalent-on-c-which-is-the-best-one 


set(APPLE_SILICON False)
if(APPLE)   
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(APPLE_SILICON True)
    endif()
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if (NOT APPLE_SILICON)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -fno-math-errno -mavx -mfma") 
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -fno-math-errno") 
endif()



find_package(GSL REQUIRED)
find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY ENG_LIBRARY MEX_COMPILER MAT_LIBRARY)


message(STATUS "Matlab_INCLUDE_DIRS: ${Matlab_INCLUDE_DIRS}")


include_directories(${Matlab_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

# Target names
set(TARGET_MEX_LIB mex)



# Libraries
set(LIBS
    GSL::gsl
)



# Source files
set(DIR_MEX_SRC src)


message(STATUS "CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")

file(GLOB MEX_FILES CONFIGURE_DEPENDS ${DIR_MEX_SRC}/*.cpp ${DIR_MEX_SRC}/*.c)


message(STATUS "LIBS=${LIBS}")
message(STATUS "TARGET_LIB=${TARGET_LIB}")

message(STATUS "Matlab_ROOT_DIR=${Matlab_ROOT_DIR}")

foreach(mexfilepath IN LISTS MEX_FILES)

    get_filename_component(mexfilename ${mexfilepath} NAME_WE)

    set(TARGET_MEX "${mexfilename}")
    message(STATUS "TARGET_MEX=${TARGET_MEX}")
    file(GLOB_RECURSE MEX_SOURCES 
        ${mexfilepath} 
        "${Matlab_ROOT_DIR}/extern/version/cpp_mexapi_version.cpp")
    add_library(${TARGET_MEX} SHARED ${MEX_SOURCES})
    target_include_directories(${TARGET_MEX} PUBLIC
        ${Matlab_INCLUDE_DIRS}
        src 
        mex )
    
    target_link_libraries(${TARGET_MEX} PUBLIC ${LIBS} ${TARGET_LIB} ${Matlab_LIBRARIES})
    target_compile_options(${TARGET_MEX} PRIVATE -Wall -Wextra -Wno-sign-compare)

    set_target_properties(${TARGET_MEX} PROPERTIES SUFFIX ".${Matlab_MEX_EXTENSION}")
    set_target_properties(${TARGET_MEX} PROPERTIES PREFIX "")
    
endforeach()

