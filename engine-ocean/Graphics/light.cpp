#include "light.h"


Light::Light(LightType type, glm::vec3 lightData, glm::vec3 lightColor):
    m_lightType(type),
    m_lightColor(lightColor),
    m_lightFunction(glm::vec3(1))
{
    switch(type){
        case LightType::POINT:
            m_lightPos = lightData;
            break;
        case LightType::DIRECTIONAL:
            m_lightDir = lightData;
            break;
    }
}

Light::~Light(){
    
}

void Light::setPos(glm::vec3 newPos){
    m_lightPos = newPos;
}

glm::vec3 Light::getPos(){
    return m_lightPos;
}

void Light::translate(glm::vec3 delta){
    m_lightPos += delta;
}

void Light::setAttenuation(glm::vec3 attenuation){
    m_lightFunction = attenuation;
}

glm::vec3 Light::getAttenuation(){
    return m_lightFunction;
}

void Light::setDir(glm::vec3 newDir){
    m_lightDir = newDir;
}

glm::vec3 Light::getDir(){
    return m_lightDir;
}

void Light::rotate(float angle, glm::vec3 axis){
    m_lightDir = glm::vec3(glm::rotate(glm::mat4(1), angle, axis)*glm::vec4(m_lightDir, 0));
}

void Light::setColor(glm::vec3 newColor){
    m_lightColor = newColor;
}

glm::vec3 Light::getColor(){
    return m_lightColor;
}

void Light::setType(LightType newType){
    m_lightType = newType;
}

LightType Light::getType(){
    return m_lightType;
}
