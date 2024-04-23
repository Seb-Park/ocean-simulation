#include "Application.h"
#include "Game/menuscreen.h"

Application::Application()
{
    game = std::make_shared<GameplayScreen>(m_input_map);//new GameplayScreen(m_input_map);
    menu = std::make_shared<MenuScreen>(m_input_map);
    activeScreen = game;

    createKeyInput(GLFW_KEY_M);
    createKeyInput(GLFW_KEY_B);
    createKeyInput(GLFW_KEY_R);

}



void Application::createKeyInput(int inputVal){
    Input input;
    input.inputVal = inputVal;
    m_input_map[inputVal] = input;
}

void Application::inactivateOtherKeys(int keyVal){
    for (auto &keypair : m_input_map){
        if (keypair.second.inputVal != keyVal){
            m_input_map.at(keypair.first).isActive = false;
        }
    }
}



void Application::update(double deltaTime){
    activeScreen->update(deltaTime);
}

void Application::draw(){
    if (m_input_map.at(GLFW_KEY_M).isActive){
        activeScreen = menu;
    } else {
        activeScreen = game;
    }

    activeScreen->draw();
}

void Application::keyEvent(int key, int action){
    switch(key){
        case GLFW_KEY_M:
            if (action == GLFW_PRESS){
                inactivateOtherKeys(key);
                m_input_map.at(key).isActive = true;
            }
        break;
        case GLFW_KEY_B:
            if (action == GLFW_PRESS){
                inactivateOtherKeys(key);
                m_input_map.at(key).isActive = true;
            }
        break;
        case GLFW_KEY_R:
            if (action == GLFW_PRESS){
                inactivateOtherKeys(key);
                m_input_map.at(key).isActive = true;
                game = std::make_shared<GameplayScreen>(m_input_map);//new GameplayScreen(m_input_map);
                menu = std::make_shared<MenuScreen>(m_input_map);
                activeScreen = game;
            }
        break;
        default:
        break;
    }

    activeScreen->keyEvent(key, action);

}

void Application::mousePosEvent(double xpos, double ypos){

    activeScreen->mousePosEvent(xpos, ypos);

}

void Application::mouseButtonEvent(int button, int action){
   activeScreen->mouseButtonEvent(button, action);

}

void Application::scrollEvent(double distance){
    activeScreen->scrollEvent(distance);

}

void Application::framebufferResizeEvent(int width, int height){
    activeScreen->framebufferResizeEvent(width, height);
}

void Application::windowResizeEvent(int width, int height){
    activeScreen->windowResizeEvent(width, height);
}
