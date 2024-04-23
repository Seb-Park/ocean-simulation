#ifndef AISYSTEM_H
#define AISYSTEM_H


// where is player?
#include "Game/GameObjects/GameObject.h"
#include <map>
#include <memory>
#include <glm/glm.hpp>


struct PosData{
    glm::vec3 currPos;
    glm::vec3 setToPos;
};

// did player jump? did player slash?
struct ConditionData{
    bool conditionTrue = false;
};

// BLACKBOARD
struct BlackboardData{
    PosData locationData;
    std::map<std::string, ConditionData> conditionData;
    std::shared_ptr<GameObject> environment;
};


class AISystem
{
public:
    AISystem();
};

#endif // AISYSTEM_H
