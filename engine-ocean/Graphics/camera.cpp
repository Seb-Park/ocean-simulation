#include "camera.h"
#include <iostream>

Camera::Camera(int width, int height, glm::vec3 pos, glm::vec3 look, 
    glm::vec3 up, float fov, float nearPlane, float farPlane):
    m_width(width),
    m_height(height),
    m_pos(pos),
    m_look(look),
    m_up(up),
    m_fov(fov),
    m_aspect(width/float(height)),
    m_near(nearPlane),
    m_far(farPlane)
{
    calculateProjection();
    calculateView();
}

Camera::~Camera(){
    
}

glm::mat4 Camera::getProjection(){
    return m_proj;
}

glm::mat4 Camera::getView(){
    return m_view;
}

void Camera::resize(int width, int height){
    m_width = width;
    m_height = height;
    m_aspect = float(width)/float(height);

    calculateProjection();
}

void Camera::translate(glm::vec3 move){
    m_pos += move;

    calculateView();
}

void Camera::setPos(glm::vec3 newPos){
    m_pos = newPos;

    calculateView();
}

glm::vec3 Camera::getPos(){
    return m_pos;
}

void Camera::rotate(float angle, glm::vec3 axis){
    glm::mat4 lookRotation = glm::rotate(glm::mat4(1), angle, axis);
    glm::vec3 tempLook = glm::vec3(lookRotation * glm::vec4(m_look, 0));
    if(glm::cross(tempLook, m_up) != glm::vec3(0)){
        m_look = tempLook;
       calculateView();
    }
}

void Camera::setLook(glm::vec3 newLook){
//    std::cout << "up: " << m_up.x << "," << m_up.y << "," << m_up.z << std::endl;
//    std::cout << "look: " << newLook.x << "," << newLook.y << "," << newLook.z << std::endl;
    if(glm::cross(newLook, m_up) == glm::vec3(0)){
        std::cout<<"Error: Look vector cannot be parallel to up vector!"<<std::endl;
        return;
    }
    m_look = newLook;

    calculateView();
}

glm::vec3 Camera::getLook(){
    return m_look;
}

void Camera::setUp(glm::vec3 newUp){
    if(glm::cross(newUp, m_look) == glm::vec3(0)){
        std::cout<<"Error: Up vector cannot be parallel to look vector!"<<std::endl;
        return;
    }
    m_up = newUp;

    calculateView();
}

glm::vec3 Camera::getUp(){
    return m_up;
}

void Camera::calculateProjection(){
    m_proj = glm::perspective(m_fov, m_aspect, m_near, m_far);
}

void Camera::calculateView(){
    // eye, lookat, up
    m_view = glm::lookAt(m_pos, m_look, m_up);
}

//

glm::vec3 Camera::getViewDirection(){
    return -glm::transpose(m_view)[2];
}

glm::vec3 Camera::getRight(){
    return glm::transpose(m_view)[0];
}

glm::vec3 Camera::getLookAt(){
    return m_look;
}

int Camera::getWidth(){
    return m_width;
}

int Camera::getHeight(){
    return m_height;
}
