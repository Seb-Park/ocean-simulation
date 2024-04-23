#include "transformcomponent.h"

TransformComponent::TransformComponent(std::shared_ptr<ModelTransform> mt,
                                       std::string entity_id,
                                       std::map<std::string,
                                       BlackboardData>& global_blackboard, bool isAI):
    m_model_transform(mt),
    m_global_blackboard(global_blackboard),
    isAIObject(isAI)
{
    old_pos = mt->getPos();
    estimated_final_pos = mt->getPos();
    m_entity_id = entity_id;

    // initialize blackboard data for entity
    BlackboardData data;
    m_global_blackboard.insert(std::pair<std::string, BlackboardData>(m_entity_id, data));
    m_global_blackboard[m_entity_id].locationData.currPos = old_pos;
}

std::shared_ptr<ModelTransform> TransformComponent::getMT(){
    return m_model_transform;
}
std::vector<std::shared_ptr<ModelTransform>> TransformComponent::getAllMT(){
     return m_all_model_transforms;
}

bool TransformComponent::hasMultipleMT(){
    return multipleMT;
}

void TransformComponent::translate(const glm::vec3 &delta){
    m_model_transform->translate(delta);

    // update old pos for collisions

    old_pos = m_model_transform->getPos();

    if (isAIObject){
        m_global_blackboard[m_entity_id].locationData.currPos = old_pos;
    }
}

void TransformComponent::setPos(const glm::vec3 &new_pos){
    m_model_transform->setPos(new_pos);

    // update old pos for collisions
    old_pos = new_pos;

    if (isAIObject){
        m_global_blackboard[m_entity_id].locationData.currPos = old_pos;
    }
}

glm::vec3 TransformComponent::getPos(){
    return m_model_transform->getPos();
}

void TransformComponent::setScale(const glm::vec3 &scale){
    m_model_transform->setScale(scale);
}

glm::vec3 TransformComponent::getScale(){
    return m_model_transform->getScale();
}

float TransformComponent::getYRotationAngle(){
    glm::mat4 rotMat = m_model_transform->getRotation();
    float ry0 = rotMat[0][2];
//    float ry1 = rotMat[1][2];
//    float ry2 = rotMat[2][2];

//    float sign = 1.f;
//    if (ry1*ry2 < 0) sign = -1.f;

//    float angley = std::atan2(-ry0, sign*pow((pow(ry1,2) + pow(ry2,2)), .5));
    float angley = -std::asin(ry0);
    float angley2 = M_PI - angley;

    //std::cout << "angle 1: " << angley << std::endl;
    std::cout << "angle 2: " << angley2 << std::endl;

    return angley;
}

void TransformComponent::setRotation(float angle, glm::vec3 axis){
    m_model_transform->rotate(angle, axis);
}




