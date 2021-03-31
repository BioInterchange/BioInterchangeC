cmake_minimum_required(VERSION 3.0.2)

set(VERSION_NUMBER 2.0.5)
set(BUILD_NUMBER 139)
set(BIOINTERCHANGE_VERSION_CMAKE ${VERSION_NUMBER}.${BUILD_NUMBER})
set(BIOINTERCHANGE_VERSION_CMPLD ${VERSION_NUMBER}+${BUILD_NUMBER})
add_definitions(-DBIOINTERCHANGE_VERSION="${BIOINTERCHANGE_VERSION_CMPLD}")

project(biointerchange VERSION ${BIOINTERCHANGE_VERSION_CMAKE} LANGUAGES C CXX)

find_program(BASH bash HINTS /bin)
find_program(BREW brew HINTS /usr/bin /usr/local/bin /home/linuxbrew/.linuxbrew/bin)

execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions curl | sed -E 's/.* //'" OUTPUT_VARIABLE CURL_VER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions googletest | sed -E 's/.* //'" OUTPUT_VARIABLE GTEST_VER OUTPUT_STRIP_TRAILING_WHITESPACE)
if(APPLE)
    execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions libiconv | sed -E 's/.* //'" OUTPUT_VARIABLE ICONV_VER OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(APPLE)
execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions openssl@1.1 | sed -E 's/.* //' | sed -E 's/\\.[[:digit:]]+[[:alnum:]]?$//'" OUTPUT_VARIABLE OPENSSL_VER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions openssl@1.1 | sed -E 's/.* //' | sed -E 's/[[:digit:]]+\\.[[:digit:]]+\\.//'" OUTPUT_VARIABLE OPENSSL_PATCH OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions python@3.9 | sed -E 's/.* //' | sed -E 's/\\.[[:digit:]]+_?[[:alnum:]]*$//'" OUTPUT_VARIABLE PYTHON_VER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions python@3.9 | sed -E 's/.* //' | sed -E 's/(^[[:digit:]]+\\.[[:digit:]]\\.|_[[:alnum:]]*$)//g'" OUTPUT_VARIABLE PYTHON_PATCH OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${BASH} -l -c "${BREW} list --versions python@3.9 | sed -E 's/.* //' | sed -E 's/^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+//'" OUTPUT_VARIABLE PYTHON_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

# Prevent `make install` to re-build all external projects yet again.
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif(APPLE)

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

# NOTE: The brew build will not see LibDocument's CPython build and other build
#       will not see `src/python` since the external project wull not be included.
include_directories("${CMAKE_SOURCE_DIR}/include"
    "/usr/local/Cellar/curl/${CURL_VER}/include"
    "/usr/local/Cellar/googletest/${GTEST_VER}/include/gtest"
    "/usr/local/Cellar/openssl@${OPENSSL_VER}/${OPENSSL_VER}.${OPENSSL_PATCH}/include"
    "/usr/local/Cellar/python@${PYTHON_VER}/${PYTHON_VER}.${PYTHON_PATCH}${PYTHON_SUFFIX}/Frameworks/Python.framework/Versions/${PYTHON_VER}/include/python${PYTHON_VER}"
    "${CMAKE_SOURCE_DIR}/src/python"
    "${CMAKE_SOURCE_DIR}/src/python/Include"
    "${CMAKE_SOURCE_DIR}/src/libdocument/include"
    "${CMAKE_SOURCE_DIR}/src/libdocument/src/python"
    "${CMAKE_SOURCE_DIR}/src/libdocument/src/python/Include")

link_directories(lib)

add_library(libdocument-static STATIC IMPORTED)
set_target_properties(libdocument-static PROPERTIES
    IMPORTED_LOCATION src/libdocument/lib/libdocument-static${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(biointerchange-static STATIC ${biointerchange_SRC})
set_target_properties(biointerchange-static PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY lib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE lib
    IMPORT_LOCATION lib/${CMAKE_CFG_INTDIR}/biointerchange${CMAKE_STATIC_LIBRARY_SUFFIX})

find_library(CURL NAMES libcurl${CMAKE_STATIC_LIBRARY_SUFFIX}
    HINTS "/usr/local/Cellar/curl/${CURL_VER}/lib" "/usr/local/lib" "/usr/lib")

# NOTE: Windows builds are not supported as such, but these names
#       can only be used with OpenSSL 1.1 and later.
find_library(CRYPTO_LIB NAMES libcrypto${CMAKE_STATIC_LIBRARY_SUFFIX}
    HINTS "/usr/local/Cellar/openssl@${OPENSSL_VER}/${OPENSSL_VER}.${OPENSSL_PATCH}/lib" "/usr/local/lib" "/usr/lib")
find_library(SSL_LIB NAMES libssl${CMAKE_STATIC_LIBRARY_SUFFIX}
    HINTS "/usr/local/Cellar/openssl@${OPENSSL_VER}/${OPENSSL_VER}.${OPENSSL_PATCH}/lib" "/usr/local/lib" "/usr/lib")

add_dependencies(biointerchange-static libdocument-static)
if(APPLE)
    find_library(PYTHON_LIB NAMES libpython${PYTHON_VER}${CMAKE_STATIC_LIBRARY_SUFFIX}
        HINTS "/usr/local/Cellar/python@${PYTHON_VER}/${PYTHON_VER}.${PYTHON_PATCH}${PYTHON_SUFFIX}/Frameworks/Python.framework/Versions/${PYTHON_VER}/lib/python${PYTHON_VER}/config-${PYTHON_VER}-darwin" "/usr/local/lib" "/usr/lib")
    find_library(CORE_FOUNDATION_LIB CoreFoundation)
    find_library(ICONV_LIB NAMES libiconv${CMAKE_STATIC_LIBRARY_SUFFIX}
        HINTS "/usr/local/Cellar/libiconv/${ICONV_VER}/lib" "/usr/local/lib")
    target_link_libraries(biointerchange-static libdocument-static ${CORE_FOUNDATION_LIB} ${ICONV_LIB} ${SSL_LIB} ${CRYPTO_LIB} ${CURL_LIB} ${PYTHON_LIB} ${INTL_LIB})
else(APPLE)
    find_library(PYTHON_LIB NAMES libpython${PYTHON_VER}${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "/usr/local/Cellar/python@${PYTHON_VER}/${PYTHON_VER}.${PYTHON_PATCH}${PYTHON_SUFFIX}/Frameworks/Python.framework/Versions/${PYTHON_VER}/lib/python${PYTHON_VER}/config-${PYTHON_VER}-darwin" "/usr/local/lib" "/usr/lib")
    target_link_libraries(biointerchange-static libdocument-static pthread ${SSL_LIB} ${CRYPTO_LIB} ${CURL_LIB} ${PYTHON_LIB} dl m util)
endif(APPLE)

add_executable(biointerchange src/biointerchange/main.c)
add_dependencies(biointerchange biointerchange-static)

add_executable(tests ${test_SRC})
add_dependencies(tests biointerchange-static)

target_link_libraries(biointerchange biointerchange-static)

install(PROGRAMS biointerchange TYPE BIN)

# As suggested in the googletest README: link libgtest statically:
find_library(SSL_LIB NAMES libssl${CMAKE_STATIC_LIBRARY_SUFFIX}
    HINTS "/usr/local/Cellar/googletest/${GTEST_VER}/lib" "/usr/local/lib" "/usr/lib")
target_link_libraries(tests biointerchange-static ${GTEST_LIB})

