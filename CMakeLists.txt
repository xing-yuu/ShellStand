cmake_minimum_required(VERSION 3.2)
project(ARAP_deformation)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_CXX_STANDARD 11)

find_package(Qt5
    REQUIRED
    COMPONENTS
    Core
    Gui
    Widgets
    )


option(BUILD_SHARED_LIBS "" OFF)
option(USE_SYSTEM_VTK "" ON)
option(USE_SYSTEM_ITK "" ON)

set(MOSEK_INCLUDE_DIR "C:/Program Files/Mosek/10.1/tools/platform/win64x86/h")
set(MOSEK_LIB_DIR "C:/Program Files/Mosek/10.1/tools/platform/win64x86/bin")
set(YAML_INCLUDE_DIR "D:/library/yaml/yaml-cpp/include")
set(YAML_LIB_DIR "D:/library/yaml/yaml-cpp/build")

FIND_LIBRARY(YAML_LIBRARY_RELEASE yaml-cpp.lib ${YAML_LIB_DIR}/Release)
FIND_LIBRARY(YAML_LIBRARY_DEBUG yaml-cppd.lib ${YAML_LIB_DIR}/Debug)
# FIND_LIBRARY(SequentialLineSearch_LIBRARY_debug SequentialLineSearch.lib ${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/lib/Debug)

include_directories(${YAML_INCLUDE_DIR})

link_directories(${MOSEK_LIB_DIR})

set(SOURCES
    main.cpp
    maths_funcs.cpp
    universal.cpp
    ThickeningRange.cpp
    Evaluation.cpp
    SelfIntersect.cpp
    Deformation.cpp
    Visibility.cpp
    Saliency.cpp
    MOSEKOptimize.cpp
    ShellThickening.cpp
    BayesianOptimization.cpp
    ModelOperation.cpp
)

# file(GLOB
# 	headers
	# ${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/include/sequential-line-search/*.hpp     //ttt
# )
file(GLOB
	sources
	${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/src/*.cpp
)

# bayesian optimization
# FIND_LIBRARY(SequentialLineSearch_LIBRARY_release SequentialLineSearch.lib ${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/lib/Release)
# FIND_LIBRARY(SequentialLineSearch_LIBRARY_debug SequentialLineSearch.lib ${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/lib/Debug)


# add_executable(${PROJECT_NAME} ${SOURCES})

# add_executable(${PROJECT_NAME} ${SOURCES} ${headers} ${sources}) ///ttt
add_executable(${PROJECT_NAME} ${SOURCES} ${sources})
target_include_directories(${PROJECT_NAME} PRIVATE  ${PROJECT_SOURCE_DIR} ${headers} ${MOSEK_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/include)
target_link_libraries ( ${PROJECT_NAME} 
                        
						debug ${YAML_LIBRARY_DEBUG}
						optimized ${YAML_LIBRARY_RELEASE} )

find_package(cinolib REQUIRED)
target_link_libraries(${PROJECT_NAME} cinolib)
find_package(VTK 8.2.0 REQUIRED)
include(${VTK_USE_FILE})
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
find_package(Eigen3 REQUIRED)
include( ${EIGEN3_USE_FILE} )



target_link_libraries(
    ${PROJECT_NAME}
    ${VTK_LIBRARIES}
    "${CGAL_LIBS}"
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    )


# add_library(${PROJECT_NAME} STATIC ${headers} ${sources})


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/external/timer)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/external/mathtoolbox)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/external/nlopt)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/external/nlopt-util)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bayesian-optimization/external/parallel-util)


target_link_libraries(${PROJECT_NAME} Eigen3::Eigen mathtoolbox nlopt nlopt-util parallel-util timer mosek64_10_1)