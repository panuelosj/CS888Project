cmake_minimum_required (VERSION 3.0)
project (CS888Project)

include_directories (include)
file (GLOB SOURCES "src/*.cpp")

add_subdirectory(Dependencies/glfw)
include_directories(Dependencies/glfw/include)

add_executable(CS888Project ${SOURCES})
target_link_libraries(CS888Project ${OPENGL_gl_LIBRARY} glfw)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

target_link_libraries(CS888Project  ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} glfw ${GLFW_LIBRARIES})
