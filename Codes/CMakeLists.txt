cmake_minimum_required (VERSION 3.5)

project ("AC24HourglassSupplementaryMaterials")

FIND_PACKAGE( OpenMP REQUIRED) 
if(OPENMP_FOUND) 
message("OPENMP FOUND") 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}") 
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}") 
set(CMAKE_SHARE_LINKER_FLAGS "${CMAKE_SHARE_LINKER_FLAGS} ${OpenMP_SHARE_LINKER_FLAGS}") 
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-mavx2)
endif()

set(CMAKE_CXX_STANDARD 14)
#.exe
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#.so
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#.a
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_executable(testAlzette mainAlzetteExpDL.cpp basics.cpp basics.h)
add_executable(testChacha mainChachaExpDL.cpp basics.cpp basics.h)
add_executable(testChaskey mainChaskeyExpDL.cpp basics.cpp basics.h)
add_executable(testForro mainForroExpDL.cpp basics.cpp basics.h)
add_executable(testSalsa mainSalsaExpDL.cpp basics.cpp basics.h)
add_executable(testSipHash mainSipHashExpCor.cpp basics.cpp basics.h)
add_executable(testSpeck32 mainSpeck32ExpDL.cpp basics.cpp basics.h)
add_executable(testSpeck48 mainSpeck48ExpDL.cpp basics.cpp basics.h)
add_executable(testSpeck64 mainSpeck64ExpDL.cpp basics.cpp basics.h)
add_executable(testSpeck96 mainSpeck96ExpDL.cpp basics.cpp basics.h)
add_executable(testSpeck128 mainSpeck128ExpDL.cpp basics.cpp basics.h)

