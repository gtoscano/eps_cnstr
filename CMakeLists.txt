cmake_minimum_required(VERSION 3.14)
project(eps_cnstr)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "-O3")

# -lboost_log
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# -DBOOST_LOG_DYN_LINK
# definitions
add_definitions(-DBOOST_ALL_NO_LIB)
add_definitions(-DBOOST_LOG_DYN_LINK=1)

# -lboost_log -lboost_thread -lboost_log_setup
find_package(Boost 1.74.0 COMPONENTS thread log log_setup REQUIRED)
find_package(PkgConfig)
find_package(Parquet REQUIRED)

# Result: -L/usr/local/lib -lpthread -lprotobuf -lsnappy -lfmt -lcurl -lloki-cpp
# pkg_check_modules(PROTOBUF REQUIRED IMPORTED_TARGET protobuf)
# pkg_check_modules(SNAPPY REQUIRED IMPORTED_TARGET snappy)
pkg_check_modules(CURLPP REQUIRED IMPORTED_TARGET curlpp)
pkg_check_modules(LIBCURL REQUIRED IMPORTED_TARGET libcurl)
pkg_check_modules(IPOPT REQUIRED IMPORTED_TARGET ipopt)

#find_package(Boost 1.74.0 COMPONENTS log_setup REQUIRED)
# git clone https://github.com/redis/hiredis.git
# cd hiredis; mkdir -p build; cd build; cmake ..; make; sudo make install
# find_package(hiredis)
# find_package(redis++)

# lcurlp
#FetchContent_Declare(curlpp
#        GIT_REPOSITORY https://github.com/jpbarrette/curlpp.git
#        GIT_TAG master
#        )

# -lfmt
#FetchContent_Declare(fmt
#  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#  GIT_TAG master
#)

#loki-cpp, dependencies: fmt, curl protobuf + snappy (optional)
#FetchContent_Declare(loki-cpp
#        GIT_REPOSITORY https://github.com/mircodezorzi/loki-cpp.git
#        GIT_TAG master
#        )
# -lhiredis
#FetchContent_Declare(hiredis
#        GIT_REPOSITORY https://github.com/redis/hiredis.git
#        GIT_TAG master
#        )
# -lhiredis
#FetchContent_Declare(redis++
#        GIT_REPOSITORY https://github.com/sewenew/redis-plus-plus
#        GIT_TAG master
#        )

      #FetchContent_Declare(re2
      #        GIT_REPOSITORY https://github.com/google/re2.git
      #        GIT_TAG main
      #        )
      #FetchContent_Declare(c-ares
      #        GIT_REPOSITORY https://github.com/c-ares/c-ares
      #        GIT_TAG main
      #        )
# -lSimpleAmqpClient
#FetchContent_Declare(SimpleAmqpClient
#        GIT_REPOSITORY https://github.com/alanxz/SimpleAmqpClient
#        GIT_TAG master
#        )
#
#FetchContent_MakeAvailable(fmt hiredis redis++ SimpleAmqpClient curlpp)

set(Parquet_DIR /usr/lib/x86_64-linux-gnu/cmake/Arrow)

# -lpthread
find_package(Threads REQUIRED)
# -larrow
find_package(Arrow CONFIG REQUIRED)
# -lparquet
find_package(Parquet CONFIG REQUIRED)
# lcrossguid -luuid
# sudo apt-get install uuid-dev
find_package(crossguid CONFIG REQUIRED)
# -laws-cpp-sdk-core
find_package(aws-cpp-sdk-core CONFIG REQUIRED)
# -laws-cpp-sdk-s3
find_package(aws-cpp-sdk-s3 CONFIG REQUIRED)
# find_package(Protobuf MODULE REQUIRED)


file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(eps_cnstr ${SRC_FILES})
target_include_directories(eps_cnstr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(eps_cnstr PUBLIC Boost::log Boost::log_setup hiredis redis++ fmt arrow_shared parquet_shared crossguid aws-cpp-sdk-core aws-cpp-sdk-s3 SimpleAmqpClient Threads::Threads PkgConfig::IPOPT PkgConfig::LIBCURL)
#PkgConfig::PROTOBUF PkgConfig::SNAPPY)

install(TARGETS eps_cnstr RUNTIME DESTINATION /usr/local/bin)

