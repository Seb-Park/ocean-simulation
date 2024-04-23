#pragma once

#include "glm/glm.hpp"
#include "glm/ext.hpp"

enum LightType{
    POINT,
    DIRECTIONAL
};

class Light{
public:
    Light(LightType type, glm::vec3 lightData, glm::vec3 lightColor = glm::vec3(1));
    ~Light();

    // Functions for point lights
    void setPos(glm::vec3 newPos);
    glm::vec3 getPos();
    void translate(glm::vec3 delta);
    void setAttenuation(glm::vec3 attenuation);
    glm::vec3 getAttenuation();

    // Functions for directional lights
    void setDir(glm::vec3 newDir);
    glm::vec3 getDir();
    void rotate(float angle, glm::vec3 axis);

    // Utility functions for light color
    void setColor(glm::vec3 newColor);
    glm::vec3 getColor();

    // Utility functions for the light type
    void setType(LightType newType);
    LightType getType();

private:
    LightType m_lightType;
    glm::vec3 m_lightPos;
    glm::vec3 m_lightDir;
    glm::vec3 m_lightColor;
    glm::vec3 m_lightFunction;
};