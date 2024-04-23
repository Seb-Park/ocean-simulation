#pragma once

#include "glm/glm.hpp"
#include "glm/ext.hpp"

class ModelTransform{
public:
    ModelTransform();
    ~ModelTransform();

    void setPos(glm::vec3 pos);
    void translate(glm::vec3 delta);
    glm::vec3 getPos();
    void setScale(glm::vec3 scale);
    void setScale(float scale);
    void scale(glm::vec3 scale);
    void scale(float scale);
    glm::vec3 getScale();
    void setRotation(float angle, glm::vec3 axis);
    void setRotation(glm::mat4 r);
    void rotate(float angle, glm::vec3 axis);
    void rotate(glm::mat4 r);
    glm::mat4 getRotation();
    glm::mat4 getModelMatrix();
    void clear();

    float getYRotationAngle();


private:
    glm::vec3 m_pos;
    glm::vec3 m_scale;
    glm::mat4 m_rotation = glm::mat4(1);
};
