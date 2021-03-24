cmake_minimum_required(VERSION 3.0.2)

execute_process(COMMAND git rev-list --all
    COMMAND wc -l
    COMMAND tr -d " \\n"
    OUTPUT_VARIABLE BUILD_NUMBER)
set(VERSION_NUMBER 2.0.5)
set(BIOINTERCHANGE_VERSION_CMAKE ${VERSION_NUMBER}.${BUILD_NUMBER})
set(BIOINTERCHANGE_VERSION_CMPLD ${VERSION_NUMBER}+${BUILD_NUMBER})
add_definitions(-DBIOINTERCHANGE_VERSION="${BIOINTERCHANGE_VERSION_CMPLD}")

project(biointerchange VERSION ${BIOINTERCHANGE_VERSION_CMAKE} LANGUAGES C CXX)

# Note: Currently using LibDocument's build!
set(PYTHON_VERSION "3.9")
set(PYTHON_PATCH "2")
set(PYTHON_SUFFIX "")

# Prevent `make install` to re-build all external projects yet again.
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

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
#    CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles"
    CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles"
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

ExternalProject_Add(googletest-lib
    GIT_REPOSITORY https://github.com/google/googletest.git
    #SVN_REVISION -r "700"
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/googletest
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/googletest
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND autoreconf -fi COMMAND ./configure
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)

if(APPLE)
ExternalProject_Add(curl-lib
    GIT_REPOSITORY https://github.com/bagder/curl.git
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    CONFIGURE_COMMAND cmake -DBUILD_SHARED_LIBS=OFF -DHTTP_ONLY=ON -DCURL_ZLIB=OFF -DOPENSSL_ROOT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src/openssl -DOPENSSL_INCLUDE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src/openssl/include -G "Unix Makefiles"
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)
else(APPLE)
ExternalProject_Add(curl-lib
    GIT_REPOSITORY https://github.com/bagder/curl.git
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/curl
    CONFIGURE_COMMAND cmake -DBUILD_SHARED_LIBS=OFF -DHTTP_ONLY=ON -DCURL_ZLIB=OFF -G "Unix Makefiles"
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)
endif(APPLE)

if(APPLE)
ExternalProject_Add(openssl-lib
    GIT_REPOSITORY https://github.com/openssl/openssl.git
    GIT_TAG OpenSSL_1_1_1j
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/openssl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/openssl
    CONFIGURE_COMMAND ./Configure no-shared darwin64-x86_64-cc
    BUILD_COMMAND make build_libs
    INSTALL_COMMAND true)
endif(APPLE)

if(BREW_BUILD)
# NOTE: For brew, the LibDocument CPython build cannot be re-used.
ExternalProject_Add(python
    GIT_REPOSITORY https://github.com/python/cpython.git
    GIT_TAG v${PYTHON_VERSION}.${PYTHON_PATCH}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/python
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/python
    CONFIGURE_COMMAND ./configure
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND true)
endif(BREW_BUILD)

# NOTE: The brew build will not see LibDocument's CPython build and other build
#       will not see `src/python` since the external project wull not be included.
include_directories("${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/src/googletest/include"
    "${CMAKE_SOURCE_DIR}/src/curl/include/curl"
    "${CMAKE_SOURCE_DIR}/src/python"
    "${CMAKE_SOURCE_DIR}/src/python/Include"
    "${CMAKE_SOURCE_DIR}/src/libdocument/include"
    "${CMAKE_SOURCE_DIR}/src/libdocument/src/python"
    "${CMAKE_SOURCE_DIR}/src/libdocument/src/python/Include"
    "${CMAKE_SOURCE_DIR}/src/openssl/include")

link_directories(lib)

add_library(libdocument-static STATIC IMPORTED)
set_target_properties(libdocument-static PROPERTIES
    IMPORTED_LOCATION src/libdocument/lib/libdocument-static${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libcurl-static STATIC IMPORTED)
set_target_properties(libcurl-static PROPERTIES
    IMPORTED_LOCATION src/curl/lib/libcurl${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libpython-static STATIC IMPORTED)
if(BREW_BUILD)
# NOTE: For brew we cannot rely on the CPython build from LibDocument
set_target_properties(libpython-static PROPERTIES
    IMPORTED_LOCATION src/python/libpython${PYTHON_VERSION}${PYTHON_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
else(BREW_BUILD)
set_target_properties(libpython-static PROPERTIES
    IMPORTED_LOCATION src/libdocument/src/python/libpython${PYTHON_VERSION}${PYTHON_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
endif(BREW_BUILD)

if(APPLE)
add_library(libcrypto-static STATIC IMPORTED)
set_target_properties(libcrypto-static PROPERTIES
    IMPORTED_LOCATION src/openssl/libcrypto${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(libssl-static STATIC IMPORTED)
set_target_properties(libssl-static PROPERTIES
    IMPORTED_LOCATION src/openssl/libssl${CMAKE_STATIC_LIBRARY_SUFFIX})
endif(APPLE)

add_library(biointerchange-static STATIC ${biointerchange_SRC})
set_target_properties(biointerchange-static PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY lib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE lib
    IMPORT_LOCATION lib/${CMAKE_CFG_INTDIR}/biointerchange${CMAKE_STATIC_LIBRARY_SUFFIX})

if(APPLE)
add_dependencies(biointerchange-static libdocument-static libssl-static libcrypto-static libcurl-static libpython-static)
find_library(INTL_LIB NAMES intl libintl HINTS "/usr/local/lib")
target_link_libraries(biointerchange-static libdocument-static libssl-static libcrypto-static libcurl-static libpython-static ${INTL_LIB})
else(APPLE)
add_dependencies(biointerchange-static libdocument-static libcurl-static libpython-static)
target_link_libraries(biointerchange-static libdocument-static pthread ssl crypto libcurl-static libpython-static dl m util)
endif(APPLE)

add_executable(biointerchange src/biointerchange/main.c)
add_dependencies(biointerchange biointerchange-static)

add_executable(tests ${test_SRC})
add_dependencies(tests biointerchange-static googletest-lib)

target_link_libraries(biointerchange biointerchange-static)

install(PROGRAMS biointerchange TYPE BIN)

# As suggested in the googletest README: link libgtest statically:
add_library(gtest STATIC IMPORTED)
set_target_properties(gtest PROPERTIES
    IMPORTED_LOCATION src/googletest/lib/.libs/libgtest${CMAKE_STATIC_LIBRARY_SUFFIX})
target_link_libraries(tests biointerchange-static gtest)

