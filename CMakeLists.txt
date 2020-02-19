cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(KFCmd5Pi VERSION 1.0)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
find_package(Eigen3 REQUIRED NO_MODULE)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(CCGO REQUIRED)
find_package(KFBase REQUIRED)
find_package(KFCmd REQUIRED)
find_package(MEUtils REQUIRED)
find_package(MElCalc REQUIRED)
include_directories(${CCGO_INCLUDE_DIRS})
include_directories(${KFBase_INCLUDE_DIRS})
include_directories(${KFCmd_INCLUDE_DIRS})
include_directories(${MEUtils_INCLUDE_DIRS})
include_directories(${MElCalc_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/env.sh.in env.sh)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(kf5pi-run ${SOURCES})
target_include_directories(kf5pi-run
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(kf5pi-run
  PUBLIC
  ${KFCmd_LIBRARIES}
  PRIVATE
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${CCGO_LIBRARIES}
  ${KFBase_LIBRARIES}
  ${MEUtils_LIBRARIES}
  ${MElCalc_LIBRARIES})
install(TARGETS kf5pi-run DESTINATION bin)
install(FILES ${PROJECT_BINARY_DIR}/env.sh DESTINATION env)
