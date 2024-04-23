#include "GameplayScreen.h"

GameplayScreen::GameplayScreen(std::map<int, Input>& input_map):
    m_input_map(input_map),
    game(GameWorld(input_map))

{
    //game = GameWorld(m_input_map);
}


void GameplayScreen::mousePosEvent(double xpos, double ypos){
    game.mousePosEvent(xpos, ypos);
}

void GameplayScreen::mouseButtonEvent(int button, int action){
    game.mouseButtonEvent(button, action);
}


void GameplayScreen::scrollEvent(double distance){
    game.scrollEvent(distance);
}


void GameplayScreen::update(double deltaTime){
    game.update(deltaTime);
}

void GameplayScreen::draw(){
    game.draw();
}

void GameplayScreen::keyEvent(int key, int action){
    game.keyEvent(key, action);
}

void GameplayScreen::framebufferResizeEvent(int width, int height){
    game.framebufferResizeEvent(width, height);
}

void GameplayScreen::windowResizeEvent(int width, int height){
   game.windowResizeEvent(width, height);
}
