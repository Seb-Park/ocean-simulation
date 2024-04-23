#include "GameWorld.h"
#include "Game/Environment/environmentsystem.h"
#include "Game/Systems/Inventory/inventorysystem.h"
#include "Game/Systems/Pathfinding/aimovementsystem.h"
#include "Game/Systems/camerasystem.h"
#include "Game/Systems/charactercontrollersystem.h"
#include "Game/Systems/CollisionSystems/collisionsystem.h"
#include "Game/Systems/drawsystem.h"
#include "Game/Systems/objectcreationsystem.h"
#include "Game/Systems/physicssystem.h"
#include "Game/Systems/UI/uisystem.h"



GameWorld::GameWorld(std::map<int, Input>& input_map) :
    camera(std::make_shared<Camera>()),
    m_input_map(input_map)
{

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    m_game_object_names = std::vector<std::string>();
    initializeInputMap();

    addSystem<ObjectCreationSystem>(std::make_shared<ObjectCreationSystem>(m_gameobjects,
                                                                           m_dynamic_gameobjects,
                                                                           m_rigid_gameobjects,
                                                                           m_global_blackboard,
                                                                           m_lootables));
    addSystem<DrawSystem>(std::make_shared<DrawSystem>(m_gameobjects));
    addSystem<PhysicsSystem>(std::make_shared<PhysicsSystem>(m_dynamic_gameobjects, m_global_blackboard, m_lootables));
  //  addSystem<CharacterControllerSystem>(std::make_shared<CharacterControllerSystem>(m_gameobjects, camera, m_input_map, m_global_blackboard));

//    addSystem<CollisionSystem>(std::make_shared<CollisionSystem>(m_gameobjects,
//                                                                 m_dynamic_gameobjects,
//                                                                 m_rigid_gameobjects,
//                                                                 m_input_map,
//                                                                 m_global_blackboard));
//    addSystem<AIMovementSystem>(std::make_shared<AIMovementSystem>(m_dynamic_gameobjects,
//                                                                   m_rigid_gameobjects));
    addSystem<CameraSystem>(std::make_shared<CameraSystem>(m_gameobjects, camera, m_input_map));


    // seperate systems
//    collisionSystem = std::make_shared<EnvironmentCollisionDetectionSystem>(m_dynamic_gameobjects,
//                                                                            m_rigid_gameobjects,
//                                                                            m_lootables,
//                                                                            m_global_blackboard);
    //environmentSystem = std::make_shared<EnvironmentSystem>(m_rigid_gameobjects, m_dynamic_gameobjects, camera);


}


void GameWorld::createKeyInput(int inputVal){
    Input input;
    input.inputVal = inputVal;
    m_input_map[inputVal] = input;
}

void GameWorld::initializeInputMap(){
    createKeyInput(GLFW_KEY_W);
    createKeyInput(GLFW_KEY_A);
    createKeyInput(GLFW_KEY_S);
    createKeyInput(GLFW_KEY_D);
    createKeyInput(GLFW_KEY_SPACE);
    createKeyInput(GLFW_KEY_P);
    //createKeyInput(GLFW_KEY_R);

    Input mousePress;
    mousePress.inputVal = GLFW_MOUSE_BUTTON_LEFT;
    m_input_map[GLFW_MOUSE_BUTTON_LEFT] = mousePress;

    Input scrollWheel;
    scrollWheel.inputVal = 1;
    m_input_map[1] = scrollWheel;
}

void GameWorld::mousePosEvent(double xpos, double ypos){
//    m_mousePos = glm::vec2(xpos, ypos);
    if (m_mouseIsHeldTime <= 1 && m_input_map[GLFW_MOUSE_BUTTON_LEFT].isActive){
         m_mouseIsHeldTime ++;
    }

    if (m_mouseIsHeldTime > 1 && m_input_map[GLFW_MOUSE_BUTTON_LEFT].isActive){
         m_input_map[GLFW_MOUSE_BUTTON_LEFT].isHeld = true;
         m_mouseIsHeldTime ++;
    }

    if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].isActive){

    }

    for (auto &system : m_systems){
        system.second->mousePosEvent(xpos, ypos);
    }


}

void GameWorld::scrollEvent(double distance){
    for (auto &system : m_systems){
        system.second->scrollEvent(distance);
    }
}

void GameWorld::update(double deltaTime){
    for (auto &system : m_systems){
        system.second->update(deltaTime);

        // detect collision after everything
       // collisionSystem->update(deltaTime);
    }

   // environmentSystem->update(deltaTime);


}

void GameWorld::draw(){

    // then render normally
    Global::graphics.setClearColor(glm::vec3(.77f, .85f, .99f));
    Global::graphics.clearScreen(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Global::graphics.bindShader("phong");
    Global::graphics.setCameraData(camera);

    Global::graphics.setGlobalData(glm::vec3(.5f));
    Global::graphics.setClipPlane(glm::vec4(0, -1, 0, 100000));

    for (auto &system : m_systems){
        system.second->draw();
    }

    //Global::graphics.setCameraData(camera);
    //environmentSystem->draw();

}


void GameWorld::mouseButtonEvent(int button, int action){
    if (button == GLFW_MOUSE_BUTTON_LEFT){
        if (action == GLFW_RELEASE){
            if ( m_input_map[button].isActive && !m_input_map[button].isHeld){
                // we know that button was released after one click (not hold), thus a click was made
                m_input_map[button].isClicked = true;
            }
                m_input_map[button].isActive = false;
                m_input_map[button].isHeld = false;
                m_mouseIsHeldTime = 0;
                m_input_map[button].checkClickTime = 0;
                return;
        }

        if (action == GLFW_PRESS &&  !m_input_map[button].isHeld && !m_input_map[button].isActive){
            m_input_map[button].isActive = true;
            return;
        }
    }
}

void GameWorld::keyEvent(int key, int action){
    switch(key){
        case GLFW_KEY_W:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
        case GLFW_KEY_A:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
        case GLFW_KEY_S:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
        case GLFW_KEY_D:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
        case GLFW_KEY_SPACE:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
        case GLFW_KEY_P:
            if (action == GLFW_PRESS){
                m_input_map.at(key).isActive = true;
            }
            if (action == GLFW_RELEASE){
                m_input_map.at(key).isActive = false;
            }
        break;
    default:
        break;
    }
}

void GameWorld::framebufferResizeEvent(int width, int height){
    Global::graphics.setFramebufferSize(glm::ivec2(width, height));
}

void GameWorld::windowResizeEvent(int width, int height){
    Global::graphics.setWindowSize(glm::ivec2(width, height));
    camera->resize(Global::graphics.getWindowSize().x, Global::graphics.getWindowSize().y);
}
