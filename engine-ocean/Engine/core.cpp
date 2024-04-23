#include "core.h"

Core::Core()
{
}

Core::~Core(){

}

void Core::update(double deltaTime){
    m_app.update(deltaTime);

}

void Core::draw(){
    m_app.draw();
}

void Core::keyEvent(int key, int action){
    m_app.keyEvent(key, action);
}

void Core::mousePosEvent(double xpos, double ypos){

    m_app.mousePosEvent(xpos,ypos);

}

void Core::mouseButtonEvent(int button, int action){
    m_app.mouseButtonEvent(button, action);

}

void Core::scrollEvent(double distance){
    m_app.scrollEvent(distance);

}

void Core::framebufferResizeEvent(int width, int height){
    m_app.framebufferResizeEvent(width, height);
}

void Core::windowResizeEvent(int width, int height){
    m_app.windowResizeEvent(width, height);
}
