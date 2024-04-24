cmake_minimum_required(VERSION 3.16)

# Sets project name
project(ocean LANGUAGES CXX C)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Sets C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specifies required Qt components
find_package(Qt6 REQUIRED COMPONENTS Core)
find_package(Qt6 REQUIRED COMPONENTS Concurrent)
find_package(Qt6 REQUIRED COMPONENTS Xml)
find_package(Qt6 REQUIRED COMPONENTS Widgets)
find_package(Qt6 REQUIRED COMPONENTS OpenGL)
find_package(Qt6 REQUIRED COMPONENTS OpenGLWidgets)
find_package(Qt6 REQUIRED COMPONENTS Gui)

# Allows you to include files from within those directories, without prefixing their filepaths
include_directories(src)
include_directories(libs)

# Specifies .cpp and .h files to be passed to the compiler
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/mainwindow.cpp
    src/arap.cpp
    src/glwidget.cpp
    src/graphics/camera.cpp
    src/graphics/graphicsdebug.cpp
    src/graphics/meshloader.cpp
    src/graphics/shader.cpp
    src/graphics/shape.cpp
#    src/graphics/oceanshape.cpp

    src/mainwindow.h
    src/arap.h
    src/glwidget.h
    src/graphics/camera.h
    src/graphics/graphicsdebug.h
    src/graphics/meshloader.h
    src/graphics/shader.h
    src/graphics/shape.h
#    src/graphics/oceanshape.h

    util/tiny_obj_loader.h
    util/unsupportedeigenthing/OpenGLSupport

    src/ocean/ocean.h
    src/ocean/ocean.cpp

    src/ocean/ocean_alt.h src/ocean/ocean_alt.cpp
)

# GLEW: this creates its library and allows you to `#include "GL/glew.h"`
add_library(StaticGLEW STATIC glew/src/glew.c)
include_directories(${PROJECT_NAME} PRIVATE glew/include)

# Specifies libraries to be linked (Qt components, glew, etc)
target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt::Concurrent
    Qt::Core
    Qt::Gui
    Qt::OpenGL
    Qt::OpenGLWidgets
    Qt::Widgets
    Qt::Xml
    StaticGLEW
)

# This allows you to `#include "Eigen/..."`
target_include_directories(${PROJECT_NAME} PRIVATE
    Eigen
)

# Specifies other files
qt6_add_resources(${PROJECT_NAME} "Resources"
    PREFIX
        "/"
    FILES
        resources/images/uv1.png
        resources/images/kitten.png
        resources/images/hello.png
        resources/images/sand_text.jpg
        resources/images/caustic_sample.jpg
        resources/images/anamorphic.jpg
        resources/images/daniel.jpg
        resources/images/sky_clouds.png
#        resources/images/ground.png
        resources/shaders/shader.frag
        resources/shaders/shader.vert
        resources/shaders/anchorPoint.vert
        resources/shaders/anchorPoint.geom
        resources/shaders/anchorPoint.frag
        resources/shaders/texture.vert
        resources/shaders/texture.frag


)

# GLEW: this provides support for Windows (including 64-bit)
if (WIN32)
  add_compile_definitions(GLEW_STATIC)
  target_link_libraries(${PROJECT_NAME} PRIVATE
    opengl32
    glu32
  )
endif()

# Set this flag to silence warnings on Windows
if (MSVC OR MSYS OR MINGW)
  set(CMAKE_CXX_FLAGS "-Wno-volatile")
endif()

# Set this flag to silence warnings on MacOS
if (APPLE)
  set(CMAKE_CXX_FLAGS "-Wno-deprecated-volatile")
endif()
