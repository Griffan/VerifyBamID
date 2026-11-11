# Downloads pinned releases of libdeflate and htslib and builds them as static
# libraries, so results and pileup speed don't depend on whichever htslib a
# binary finds at run time. htslib's own dependencies (libcurl, libcrypto,
# zlib, bzip2, liblzma) stay dynamically linked from the system.
#
# Sets HTSLIB_VERSION, LIBDEFLATE_VERSION, HTSLIB_INCLUDE_DIRS,
# HTSLIB_LINK_LIBRARIES and HTSLIB_BUILD_TARGET for the top-level
# CMakeLists.txt.

set(HTSLIB_VERSION 1.25)
set(LIBDEFLATE_VERSION 1.26)

# htslib's configure, and the flags passed to it, split on whitespace.
if(CMAKE_BINARY_DIR MATCHES "[ \t]")
    message(FATAL_ERROR "The bundled htslib can't be built in a directory whose path contains "
        "spaces (${CMAKE_BINARY_DIR}). Use another build directory, or link an installed "
        "htslib 1.25 or later with -DUSE_SYSTEM_HTSLIB=ON.")
endif()

include(ExternalProject)
if(POLICY CMP0135)
    # Give extracted files the extraction time, so changing a pinned URL
    # rebuilds everything that depends on it.
    cmake_policy(SET CMP0135 NEW)
endif()

# htslib's dependencies, found here so that a missing one fails at configure
# time rather than part-way through the build. Imported targets are used
# where CMake 3.10 provides them. The top-level CMakeLists.txt finds zlib.
set(missing_dependencies "")
# curl's own CMake package can describe libcurl as a target rather than as
# files, which htslib's configure couldn't use.
set(CURL_NO_CURL_CMAKE ON)
find_package(CURL)
if(NOT CURL_FOUND)
    list(APPEND missing_dependencies libcurl)
endif()
# macOS signs S3 requests with CommonCrypto from libSystem; elsewhere htslib
# needs libcrypto.
if(NOT APPLE)
    find_package(OpenSSL)
    if(NOT OPENSSL_FOUND)
        list(APPEND missing_dependencies "OpenSSL (libcrypto)")
    endif()
endif()
find_package(BZip2)
if(NOT BZIP2_FOUND)
    list(APPEND missing_dependencies bzip2)
endif()
# htslib needs only the library: without lzma.h, which the macOS SDK lacks,
# it compiles against its own declarations of the liblzma API.
find_library(LIBLZMA_LIBRARY lzma)
if(NOT LIBLZMA_LIBRARY)
    list(APPEND missing_dependencies "liblzma (xz)")
endif()
find_package(Threads)
if(NOT Threads_FOUND)
    list(APPEND missing_dependencies "POSIX threads")
endif()
if(missing_dependencies)
    string(REPLACE ";" ", " missing_dependencies "${missing_dependencies}")
    message(FATAL_ERROR "The bundled htslib needs development files for: ${missing_dependencies}. "
        "Install them (see README.md), or link an installed htslib 1.25 or later with "
        "-DUSE_SYSTEM_HTSLIB=ON.")
endif()

# The sub-builds must produce archives that link into this build, so they get
# its C and linker flags and, on macOS, its architecture, deployment target
# and SDK.
get_property(user_c_flags CACHE CMAKE_C_FLAGS PROPERTY VALUE)
get_property(user_linker_flags CACHE CMAKE_EXE_LINKER_FLAGS PROPERTY VALUE)
set(libdeflate_platform_args "")
set(htslib_platform_flags "")
if(APPLE)
    list(LENGTH CMAKE_OSX_ARCHITECTURES architecture_count)
    if(architecture_count GREATER 1)
        message(FATAL_ERROR "The bundled htslib can't be built for several architectures at once, "
            "because its configure checks for one CPU's instructions. Build each architecture "
            "separately, or link an installed universal htslib with -DUSE_SYSTEM_HTSLIB=ON.")
    endif()
    foreach(setting CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
        if(${setting})
            list(APPEND libdeflate_platform_args "-D${setting}=${${setting}}")
        endif()
    endforeach()
    if(CMAKE_OSX_ARCHITECTURES)
        string(APPEND htslib_platform_flags " -arch ${CMAKE_OSX_ARCHITECTURES}")
    endif()
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
        string(APPEND htslib_platform_flags " -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    endif()
    if(IS_ABSOLUTE "${CMAKE_OSX_SYSROOT}")
        string(APPEND htslib_platform_flags " -isysroot ${CMAKE_OSX_SYSROOT}")
    endif()
endif()

# Each version installs to its own prefix, so changing a version changes every
# include path and link line that uses it, and everything is rebuilt.
set(LIBDEFLATE_PREFIX "${CMAKE_BINARY_DIR}/bundled/libdeflate-${LIBDEFLATE_VERSION}")
set(HTSLIB_PREFIX "${CMAKE_BINARY_DIR}/bundled/htslib-${HTSLIB_VERSION}")
set(LIBDEFLATE_LIBRARY "${LIBDEFLATE_PREFIX}/lib/libdeflate.a")
set(HTSLIB_LIBRARY "${HTSLIB_PREFIX}/lib/libhts.a")

# Under a Makefile generator the sub-builds run through $(MAKE), which shares
# the parent make's jobserver; a make started any other way inside a parallel
# build loses it, which older GNU make treats as an error. Under other
# generators libdeflate builds with this build's own tool, and htslib with make
# and one job per core.
if(CMAKE_GENERATOR MATCHES "Make")
    set(LIBDEFLATE_BUILD_COMMAND "$(MAKE)" install)
    set(HTSLIB_MAKE "$(MAKE)")
else()
    set(LIBDEFLATE_BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install)
    find_program(GNU_MAKE_PROGRAM NAMES gmake make)
    if(NOT GNU_MAKE_PROGRAM)
        message(FATAL_ERROR "Building the bundled htslib needs make")
    endif()
    include(ProcessorCount)
    ProcessorCount(HTSLIB_JOBS)
    if(HTSLIB_JOBS EQUAL 0)
        set(HTSLIB_JOBS 1)
    endif()
    set(HTSLIB_MAKE ${GNU_MAKE_PROGRAM} -j${HTSLIB_JOBS})
endif()

# Both projects install during their build step, the step that declares the
# installed archive, so Ninja relinks in the same run when an archive changes.
ExternalProject_Add(libdeflate_external
    URL https://github.com/ebiggers/libdeflate/releases/download/v${LIBDEFLATE_VERSION}/libdeflate-${LIBDEFLATE_VERSION}.tar.gz
    URL_HASH SHA256=125856d4656e0feab660f94842f835923410c9281fedbcee64598c918da42b5a
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${LIBDEFLATE_PREFIX}
        # GNUInstallDirs would pick lib64 on some Linux distributions.
        -DCMAKE_INSTALL_LIBDIR=lib
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        "-DCMAKE_C_FLAGS=${user_c_flags}"
        # The sub-build uses this build's generator, so it needs the same tool.
        -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DLIBDEFLATE_BUILD_SHARED_LIB=OFF
        -DLIBDEFLATE_BUILD_GZIP=OFF
        ${libdeflate_platform_args}
    BUILD_COMMAND ${LIBDEFLATE_BUILD_COMMAND}
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${LIBDEFLATE_LIBRARY})

# htslib's configure looks for the same dependencies CMake found above, after
# the libdeflate prefix so that a libdeflate in the same place can't shadow
# ours.
set(HTSLIB_CPPFLAGS "-I${LIBDEFLATE_PREFIX}/include")
set(HTSLIB_LDFLAGS "-L${LIBDEFLATE_PREFIX}/lib")
set(dependency_include_dirs ${CURL_INCLUDE_DIRS} ${OPENSSL_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS}
    ${BZIP2_INCLUDE_DIR})
if(dependency_include_dirs)
    list(REMOVE_DUPLICATES dependency_include_dirs)
endif()
foreach(dir ${dependency_include_dirs})
    list(FIND CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "${dir}" implicit)
    if(implicit EQUAL -1)
        set(HTSLIB_CPPFLAGS "${HTSLIB_CPPFLAGS} -I${dir}")
    endif()
endforeach()
set(dependency_library_dirs "")
foreach(library ${CURL_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${ZLIB_LIBRARIES}
        ${BZIP2_LIBRARIES} ${LIBLZMA_LIBRARY})
    if(IS_ABSOLUTE "${library}")
        get_filename_component(dir "${library}" DIRECTORY)
        list(APPEND dependency_library_dirs "${dir}")
    endif()
endforeach()
if(dependency_library_dirs)
    list(REMOVE_DUPLICATES dependency_library_dirs)
endif()
foreach(dir ${dependency_library_dirs})
    list(FIND CMAKE_C_IMPLICIT_LINK_DIRECTORIES "${dir}" implicit)
    if(implicit EQUAL -1)
        set(HTSLIB_LDFLAGS "${HTSLIB_LDFLAGS} -L${dir}")
    endif()
endforeach()
set(HTSLIB_CPPFLAGS "${HTSLIB_CPPFLAGS} $ENV{CPPFLAGS}")
set(HTSLIB_LDFLAGS "${HTSLIB_LDFLAGS}${htslib_platform_flags} ${user_linker_flags}")

# Every optional feature that configure would otherwise detect is pinned, so
# libhts.a always needs the same libraries, and a missing dependency fails
# configure instead of quietly dropping remote file access or a CRAM codec.
# Plugins stay off so libcurl, S3 and GCS support are compiled into
# libhts.a. These flags and HTSLIB_LINK_LIBRARIES below must stay in step.
# Only the static library and headers are built and copied; the shared
# library and htslib's programs aren't needed.
ExternalProject_Add(htslib_external
    DEPENDS libdeflate_external
    URL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
    URL_HASH SHA256=158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --prefix=${HTSLIB_PREFIX}
        --enable-libcurl
        --enable-s3
        --enable-gcs
        --enable-lzma
        --enable-bz2
        --with-libdeflate
        --disable-plugins
        --without-external-htscodecs
        --disable-ref-cache
        CC=${CMAKE_C_COMPILER}
        "CFLAGS=-O2 -fPIC${htslib_platform_flags} ${user_c_flags}"
        "CPPFLAGS=${HTSLIB_CPPFLAGS}"
        "LDFLAGS=${HTSLIB_LDFLAGS}"
        # Copying the headers here, once per configure, keeps a rebuild of
        # libhts.a from making everything that includes them recompile.
        COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/htslib ${HTSLIB_PREFIX}/include/htslib
    BUILD_COMMAND ${HTSLIB_MAKE} lib-static
        COMMAND ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/libhts.a ${HTSLIB_LIBRARY}
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${HTSLIB_LIBRARY})

# libhts.a and libdeflate.a are linked by full path because the macOS linker
# prefers a .dylib over a .a for -l. The rest is what the configure flags
# above make libhts.a depend on.
set(HTSLIB_INCLUDE_DIRS "${HTSLIB_PREFIX}/include")
set(HTSLIB_LINK_LIBRARIES ${HTSLIB_LIBRARY} ${LIBDEFLATE_LIBRARY} ${CURL_LIBRARIES})
if(NOT APPLE)
    list(APPEND HTSLIB_LINK_LIBRARIES OpenSSL::Crypto)
endif()
list(APPEND HTSLIB_LINK_LIBRARIES ZLIB::ZLIB ${BZIP2_LIBRARIES} ${LIBLZMA_LIBRARY} Threads::Threads m)
set(HTSLIB_BUILD_TARGET htslib_external)
