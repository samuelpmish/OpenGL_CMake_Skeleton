# Projects Settings
cmake_minimum_required (VERSION 3.1)
project (opengl-skeleton)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE cpp_sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_library(opengl_tools STATIC
  src/Canvas.hpp
  src/Canvas.cpp
  src/Application.hpp
  src/Application.cpp
  src/MyApplication.hpp
  src/MyApplication.cpp
  src/glError.hpp
  src/glError.cpp
  src/Shader.hpp
  src/Shader.cpp
)

add_definitions(-DGLEW_STATIC)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
set(GLFW_USE_CHDIR OFF CACHE BOOL " " FORCE)
add_subdirectory(lib/glfw EXCLUDE_FROM_ALL)
add_subdirectory(lib/glew EXCLUDE_FROM_ALL)
add_subdirectory(lib/glm EXCLUDE_FROM_ALL)
add_subdirectory(lib/imgui)

target_link_libraries(opengl_tools PUBLIC glfw glm libglew_static imgui)

add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE opengl_tools)

add_executable(main2 src/main2.cpp)
target_link_libraries(main2 PRIVATE opengl_tools)