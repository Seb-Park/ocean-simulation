#include "menuscreen.h"
#include "Game/GameWorld.h"
#include "Graphics/global.h"

MenuScreen::MenuScreen(std::map<int, Input>& input_map):
    m_input_map(input_map)
{

}

MenuScreen::~MenuScreen(){

}


void MenuScreen::update(double deltaTime){


}

void MenuScreen::draw(){
    Global::graphics.setClearColor(glm::vec3(0.f));
    Global::graphics.clearScreen(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Global::graphics.bindShader("text");
    Global::graphics.drawUIText(Global::graphics.getFont("opensans"), "Game Menu", glm::ivec2(200, 200), AnchorPoint::TopLeft, Global::graphics.getFramebufferSize().x, 1.f, 0.1f, glm::vec3(1, 1, 1));
    Global::graphics.drawUIText(Global::graphics.getFont("opensans"), "Press 'B' to resume the game.", glm::ivec2(20, 70), AnchorPoint::TopLeft, Global::graphics.getFramebufferSize().x, 0.2f, 0.1f, glm::vec3(1, 1, 1));


}

void MenuScreen::keyEvent(int key, int action){

}

void MenuScreen::mousePosEvent(double xpos, double ypos){



}

void MenuScreen::mouseButtonEvent(int button, int action){


}

void MenuScreen::scrollEvent(double distance){


}

void MenuScreen::framebufferResizeEvent(int width, int height){
    Global::graphics.setFramebufferSize(glm::ivec2(width, height));

}

void MenuScreen::windowResizeEvent(int width, int height){
    Global::graphics.setWindowSize(glm::ivec2(width, height));

}
