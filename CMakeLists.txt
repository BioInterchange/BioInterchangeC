cmake_minimum_required(VERSION 3.0.2)

execute_process(COMMAND git rev-list --all
    COMMAND wc -l
    COMMAND tr -d " \\n"
    OUTPUT_VARIABLE BUILD_NUMBER)
set(VERSION_NUMBER 2.0.0)
set(BIOINTERCHANGE_VERSION_CMAKE ${VERSION_NUMBER}.${BUILD_NUMBER})
set(BIOINTERCHANGE_VERSION_CMPLD ${VERSION_NUMBER}+${BUILD_NUMBER})
add_definitions(-DBIOINTERCHANGE_VERSION="${BIOINTERCHANGE_VERSION_CMPLD}")

project(biointerchange VERSION ${BIOINTERCHANGE_VERSION_CMAKE} LANGUAGES C CXX)

set(PYTHON_VERSION "3.4")
set(PYTHON_PATCH "3")
set(PYTHON_SUFFIX "m")

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
    CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles"
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

ExternalProject_Add(curl-lib
    GIT_REPOSITORY https://github.com/bagder/curl.git
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    CONFIGURE_COMMAND cmake -DCURL_STATICLIB=ON -DHTTP_ONLY=ON -DCURL_ZLIB=OFF -G "Unix Makefiles"
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

ExternalProject_Add(openssl-lib
    GIT_REPOSITORY https://github.com/openssl/openssl.git
    GIT_TAG OpenSSL_1_0_1m
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/openssl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/openssl
    CONFIGURE_COMMAND ./Configure no-shared darwin64-x86_64-cc
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

ExternalProject_Add(python
    HG_REPOSITORY https://hg.python.org/cpython
    HG_TAG v${PYTHON_VERSION}.${PYTHON_PATCH}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/python
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/python
    CONFIGURE_COMMAND ./configure
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

include_directories(include)
include_directories(src/googletest/include)
include_directories(src/curl/include/curl)
include_directories(src/python)
include_directories(src/python/Include)
include_directories(src/libdocument/include)

link_directories(lib)

add_library(libdocument-static STATIC IMPORTED)
set_target_properties(libdocument-static PROPERTIES
    IMPORTED_LOCATION src/libdocument/lib/libdocument-static${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libcurl-static STATIC IMPORTED)
set_target_properties(libcurl-static PROPERTIES
    IMPORTED_LOCATION src/curl/lib/libcurl${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libpython-static STATIC IMPORTED)
set_target_properties(libpython-static PROPERTIES
    IMPORTED_LOCATION src/python/libpython${PYTHON_VERSION}${PYTHON_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libcrypto-static STATIC IMPORTED)
set_target_properties(libcrypto-static PROPERTIES
    IMPORTED_LOCATION src/openssl/libcrypto${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libssl-static STATIC IMPORTED)
set_target_properties(libssl-static PROPERTIES
    IMPORTED_LOCATION src/openssl/libssl${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(biointerchange-static STATIC ${biointerchange_SRC})
set_target_properties(biointerchange-static PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY lib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE lib
    IMPORT_LOCATION lib/${CMAKE_CFG_INTDIR}/biointerchange${CMAKE_STATIC_LIBRARY_SUFFIX})

add_dependencies(biointerchange-static libdocument-static libssl-static libcrypto-static libcurl-static libpython-static)

target_link_libraries(biointerchange-static libdocument-static libssl-static libcrypto-static libcurl-static libpython-static)

add_executable(biointerchange src/biointerchange/main.c)
add_dependencies(biointerchange biointerchange-static)

add_executable(tests ${test_SRC})
add_dependencies(tests biointerchange-static googletest-lib)

target_link_libraries(biointerchange biointerchange-static)

# As suggested in the googletest README: link libgtest statically:
add_library(gtest STATIC IMPORTED)
set_target_properties(gtest PROPERTIES
    IMPORTED_LOCATION src/googletest/lib/.libs/libgtest${CMAKE_STATIC_LIBRARY_SUFFIX})
target_link_libraries(tests biointerchange-static gtest)

