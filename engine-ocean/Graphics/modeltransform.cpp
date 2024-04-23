#include "modeltransform.h"
#include <iostream>

ModelTransform::ModelTransform():
    m_pos(glm::vec3(0)),
    m_scale(glm::vec3(1)),
    m_rotation(glm::mat4(1))
{

}

ModelTransform::~ModelTransform(){

}

void ModelTransform::setPos(glm::vec3 pos){
    m_pos = pos;
}

void ModelTransform::translate(glm::vec3 delta) {
    m_pos += delta;
}

glm::vec3 ModelTransform::getPos() {
    return m_pos;
}

void ModelTransform::setScale(glm::vec3 scale) {
    m_scale = scale;
}

void ModelTransform::setScale(float scale) {
    m_scale = glm::vec3(scale);
}

void ModelTransform::scale(glm::vec3 scale) {
    m_scale *= scale;
}

void ModelTransform::scale(float scale) {
    m_scale *= scale;
}

glm::vec3 ModelTransform::getScale() {
    return m_scale;
}

void ModelTransform::setRotation(float angle, glm::vec3 axis) {
    m_rotation = glm::rotate(glm::mat4(1), angle, axis);
}

void ModelTransform::setRotation(glm::mat4 r) {
    m_rotation = r;
}

void ModelTransform::rotate(float angle, glm::vec3 axis) {
    m_rotation = glm::rotate(m_rotation, angle, axis);
}

void ModelTransform::rotate(glm::mat4 r) {
    m_rotation = r * m_rotation;
}

glm::mat4 ModelTransform::getRotation() {
    return m_rotation;
}

//float ModelTransform::getYRotationAngle(){
//    glm::mat4 rotMat = getRotation();
//    float ry0 = rotMat[0][2];
//    float ry1 = rotMat[1][2];
//    float ry2 = rotMat[2][2];

//    float angley = std::atan2(-ry0, pow((pow(ry1, 2) + pow(ry2,2)), .5));
//    return angley;
//}

glm::mat4 ModelTransform::getModelMatrix() {
    glm::mat4 model = glm::mat4(1);
    model = glm::translate(model, m_pos);
    model = model * m_rotation;
    model = glm::scale(model, m_scale);
    return model;
}

void ModelTransform::clear(){
    m_pos = glm::vec3(0, 0, 0);
    m_scale = glm::vec3(1, 1, 1);
    m_rotation = glm::mat4(1);
}
