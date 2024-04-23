#ifndef TRANSFORMCOMPONENT_H
#define TRANSFORMCOMPONENT_H

#include "Game/Systems/aisystem.h"
#include "Graphics/global.h"
#include <memory>
#include "Component.h"

class TransformComponent : public Component
{
public:
    TransformComponent(std::shared_ptr<ModelTransform> mt,
                       std::string entity_id,
                       std::map<std::string, BlackboardData>& m_global_blackboard,
                       bool isAI = false);

    void translate(const glm::vec3 &delta);
    void setPos(const glm::vec3 &new_pos);
    void setScale(const glm::vec3 &scale);
    std::shared_ptr<ModelTransform> getMT();
    std::vector<std::shared_ptr<ModelTransform>> getAllMT();
    glm::vec3 getScale();
    //glm::vec3 getRotation();
    void setRotation(float angle, glm::vec3 axis);
    float getYRotationAngle();




    bool hasMultipleMT();
    glm::vec3 getPos();

    // used for collisions
    glm::vec3 old_pos;
    glm::vec3 estimated_final_pos;
    bool onGround = false;
    bool movingLaterally = false;
    float gravity = -25.f;
    float yVelocity = 0.f;



private:
    std::shared_ptr<ModelTransform> m_model_transform;
    std::vector<std::shared_ptr<ModelTransform>> m_all_model_transforms;
    bool multipleMT = false;
    std::map<std::string, BlackboardData>& m_global_blackboard;
    std::string m_entity_id;
    bool isAIObject = false;


};

#endif // TRANSFORMCOMPONENT_H
