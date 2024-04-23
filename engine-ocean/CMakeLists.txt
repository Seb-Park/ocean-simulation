#Sets minimum cmake version 
cmake_minimum_required(VERSION 3.14)

#Sets project name to CS1950UEngine
project(simulation-ocean)

find_package(OpenGL REQUIRED)

#Tells compiler to use c++ 20
set(CMAKE_CXX_STANDARD 20)
#Makes c++ 20 required
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Setting path macros
set(GLFW_SOURCE_DIR "External/glfw-3.3.8")
set(GLEW_SOURCE_DIR "External/glew")
set(GLM_SOURCE_DIR "External/glm-master")
set(STB_SOURCE_DIR "External/stb")
set(FREETYPE_SOURCE_DIR "External/freetype")

#GLFW settings
set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "")
set(GLFW_BUILD_TESTS OFF CACHE BOOL "")
set(GLFW_BUILD_DOCS OFF CACHE BOOL "")
set(GLFW_INSTALL OFF CACHE BOOL "")

#freetype settings
set(FT_DISABLE_ZLIB TRUE)
set(FT_DISABLE_BZIP2 TRUE)
set(FT_DISABLE_PNG TRUE)
set(FT_DISABLE_HARFBUZZ TRUE)
set(FT_DISABLE_BROTLI TRUE)


add_subdirectory(${GLFW_SOURCE_DIR})
add_subdirectory(${GLM_SOURCE_DIR})
add_subdirectory(${FREETYPE_SOURCE_DIR})

add_library(StaticGLEW STATIC ${GLEW_SOURCE_DIR}/src/glew.c)

include_directories(${GLFW_SOURCE_DIR}/include
                    ${GLFW_SOURCE_DIR}/deps
                    ${GLEW_SOURCE_DIR}/include
                    ${STB_SOURCE_DIR}
                    ${FREETYPE_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}
                    )

add_executable(${PROJECT_NAME}
    main.cpp

    Engine/core.h
    Engine/core.cpp
    Engine/window.h
    Engine/window.cpp

    Game/Application.h
    Game/Application.cpp
    Game/Screen.h
    Game/Screen.cpp
    Game/MenuScreen.h
    Game/MenuScreen.cpp
    Game/GameplayScreen.h
    Game/GameplayScreen.cpp
    Game/GameWorld.h
    Game/GameWorld.cpp

    Game/TypeMap.h


    Game/Components/Component.h
    Game/Components/Component.cpp

    Game/Components/DrawComponent.h
    Game/Components/DrawComponent.cpp
    Game/Components/TransformComponent.h
    Game/Components/TransformComponent.cpp
    Game/Components/PlayerControlComponent.h
    Game/Components/PlayerControlComponent.cpp
    Game/Components/PathfindComponent.h
    Game/Components/PathfindComponent.cpp

    Game/Components/CollisionComponents/CollisionComponent.h
    Game/Components/CollisionComponents/CollisionComponent.cpp
    Game/Components/CollisionComponents/CylinderCollider.h
    Game/Components/CollisionComponents/CylinderCollider.cpp
    Game/Components/CollisionComponents/BoundingShape.h
    Game/Components/CollisionComponents/BoundingShape.cpp

    Game/GameObjects/GameObject.h
    Game/GameObjects/GameObject.cpp

    Game/Systems/system.h
    Game/Systems/system.cpp
    Game/Systems/drawsystem.h
    Game/Systems/drawsystem.cpp
    Game/Systems/camerasystem.h
    Game/Systems/camerasystem.cpp
    Game/Systems/charactercontrollersystem.h
    Game/Systems/charactercontrollersystem.cpp
    Game/Systems/physicssystem.h
    Game/Systems/physicssystem.cpp
    Game/Systems/CollisionSystems/environmentcollisiondetectionsystem.h
    Game/Systems/CollisionSystems/environmentcollisiondetectionsystem.cpp

    Game/Systems/objectcreationsystem.h
    Game/Systems/objectcreationsystem.cpp
    Game/Systems/Inventory/inventorysystem.h
    Game/Systems/Inventory/inventorysystem.cpp
    Game/Systems/Inventory/inventoryitem.h
    Game/Systems/Inventory/inventoryitem.cpp

    Graphics/tiny_obj_loader.h

    Game/Systems/CollisionSystems/collisionsystem.h
    Game/Systems/CollisionSystems/collisionsystem.cpp
    Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h
    Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.cpp
    Game/Systems/CollisionSystems/accelerationsystem.h
    Game/Systems/CollisionSystems/accelerationsystem.cpp
    Game/Systems/CollisionSystems/BVH/bvhtree.h
    Game/Systems/CollisionSystems/BVH/bvhtree.cpp
    Game/Systems/CollisionSystems/UniformGrid/uniformgrid.h
    Game/Systems/CollisionSystems/UniformGrid/uniformgrid.cpp


    Game/Components/CollisionComponents/BoundingEllipsoid.h
    Game/Components/CollisionComponents/BoundingEllipsoid.cpp
    Game/Components/CollisionComponents/BoundingTriangle.h
    Game/Components/CollisionComponents/BoundingTriangle.cpp

    Game/Components/CollisionComponents/BoundingDynamicMesh.h
    Game/Components/CollisionComponents/BoundingDynamicMesh.cpp

    Game/Systems/Pathfinding/pathfinder.h
    Game/Systems/Pathfinding/pathfinder.cpp

    Game/Systems/Pathfinding/aimovementsystem.h
    Game/Systems/Pathfinding/aimovementsystem.cpp



    Game/Systems/AI/Actions/btaction.h
    Game/Systems/AI/Actions/btaction.cpp
    Game/Systems/AI/Conditions/btcondition.h
    Game/Systems/AI/Conditions/btcondition.cpp
    Game/Systems/AI/Actions/walkaction.h
    Game/Systems/AI/Actions/walkaction.cpp
    Game/Systems/AI/Conditions/proximitycondition.h
    Game/Systems/AI/Conditions/proximitycondition.cpp

    Game/Systems/AI/btnode.h
    Game/Systems/AI/btnode.cpp

    Game/Systems/AI/btselector.h
    Game/Systems/AI/btselector.cpp

    Game/Systems/AI/btsequence.h
    Game/Systems/AI/btsequence.cpp
    Game/Systems/aisystem.h
    Game/Systems/aisystem.cpp
    Game/Systems/AI/aibehaviorcomponent.h
    Game/Systems/AI/aibehaviorcomponent.cpp

    Game/Environment/environmentsystem.h
    Game/Environment/environmentsystem.cpp

    Game/Environment/grassenvironment.h
    Game/Environment/grassenvironment.cpp
    Game/Environment/skyboxenvironment.h
    Game/Environment/skyboxenvironment.cpp

    Game/Environment/water.h
    Game/Environment/water.cpp

    Game/Environment/Environment.h
    Game/Environment/Environment.cpp

    Game/Systems/UI/uisystem.h
    Game/Systems/UI/uisystem.cpp

    Game/Systems/UI/UITextures/uitexture.h
    Game/Systems/UI/UITextures/uitexture.cpp

    Game/Systems/UI/UITextures/uibutton.h
    Game/Systems/UI/UITextures/uibutton.cpp
    Game/Systems/UI/UITextures/uidisplay.h
    Game/Systems/UI/UITextures/uidisplay.cpp

    Game/Systems/UI/uielement.h
    Game/Systems/UI/uielement.cpp
    Game/Systems/UI/ButtonAction/buttonaction.h
    Game/Systems/UI/ButtonAction/buttonaction.cpp
    Game/Systems/UI/ButtonAction/showwindowaction.h
    Game/Systems/UI/ButtonAction/showwindowaction.cpp

    Game/Ocean/ocean.h
    Game/Ocean/ocean.cpp










    Graphics/graphics.h
    Graphics/graphics.cpp
    Graphics/global.h
    Graphics/debug.h
    Graphics/camera.h
    Graphics/camera.cpp
    Graphics/shaderloader.h
    Graphics/material.h
    Graphics/material.cpp
    Graphics/shapedata.h
    Graphics/shape.h
    Graphics/shape.cpp
    Graphics/modeltransform.h
    Graphics/modeltransform.cpp
    Graphics/light.h
    Graphics/light.cpp
    Graphics/font.h
    Graphics/font.cpp
    Graphics/textrenderer.h
    Graphics/textrenderer.cpp
    Graphics/GLWrappers/shader.h
    Graphics/GLWrappers/shader.cpp
    Graphics/GLWrappers/texture.h
    Graphics/GLWrappers/texture.cpp
    Graphics/GLWrappers/vbo.h
    Graphics/GLWrappers/vbo.cpp
    Graphics/GLWrappers/vao.h
    Graphics/GLWrappers/vao.cpp
)

target_link_libraries(${PROJECT_NAME} glfw StaticGLEW glm freetype ${OPENGL_LIBRARIES})

file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/Resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

if (WIN32)
  add_compile_definitions(GLEW_STATIC)
  target_link_libraries(${PROJECT_NAME}
      opengl32
      glu32
  )
endif()

if (UNIX AND NOT APPLE)
  target_link_libraries(${PROJECT_NAME}
      GL
  )
endif()
