cmake_minimum_required(VERSION 3.0.2)

project(biointerchange VERSION 2.0.0 LANGUAGES C CXX)

file(GLOB_RECURSE biointerchange_SRC
    FOLLOW_SYMLINKS
    "include/*.h"
    "src/biointerchange/*.c")

file(GLOB_RECURSE test_SRC
    FOLLOW_SYMLINKS
    "test/*.h"
    "test/*.cpp")

enable_testing()

include(ExternalProject)

ExternalProject_Add(libdocument-lib
    GIT_REPOSITORY https://github.com/indiedotkim/LibDocument.git
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/libdocument
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/libdocument
    CONFIGURE_COMMAND cmake -G 'Unix Makefiles'
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

ExternalProject_Add(googletest-lib
    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk
    SVN_REVISION -r "700"
    SVN_TRUST_CERT 1
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/googletest
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/googletest
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND autoreconf -fi COMMAND ./configure
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

include_directories(include)
include_directories(src/googletest/include)

link_directories(lib)

add_library(biointerchange-static STATIC ${biointerchange_SRC})
set_target_properties(biointerchange-static PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY lib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE lib
    IMPORT_LOCATION lib/${CMAKE_CFG_INTDIR}/biointerchange${CMAKE_STATIC_LIBRARY_SUFFIX})

add_executable(biointerchange src/biointerchange/main.c)

#add_executable(tests ${test_SRC})
#add_dependencies(tests biointerchange-static googletest-lib)

target_link_libraries(biointerchange biointerchange-static)

# As suggested in the googletest README: link libgtest statically:
#add_library(gtest STATIC IMPORTED)
#set_target_properties(gtest PROPERTIES
#    IMPORTED_LOCATION src/googletest/lib/.libs/libgtest${CMAKE_STATIC_LIBRARY_SUFFIX})
#target_link_libraries(tests document-static gtest)
