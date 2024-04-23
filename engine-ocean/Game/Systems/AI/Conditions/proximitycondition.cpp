#include "proximitycondition.h"
#include "Game/Components/TransformComponent.h"
#include "glm/glm.hpp"
#include <memory>

ProximityCondition::ProximityCondition(std::string entity_id,
                                       std::map<std::string, BlackboardData>& global_blackboard,
                                       float proximity):
    m_global_blackboard(global_blackboard)
{
    m_proximity = proximity;
    m_entity_id = entity_id;
    // initialize just in case
    m_global_blackboard["player"].conditionData["isJumping"].conditionTrue = false;


}

// maybe can check locations from blackboard
// pass blackboard into constructor
// struct: positiondata --> getCurrentPos, setCurrentPos
bool ProximityCondition::checkProximity(){
    // unrooted distance
    glm::vec3 aiPos = m_global_blackboard[m_entity_id].locationData.currPos;
    glm::vec3 otherPos = m_global_blackboard["player"].locationData.currPos;
    float distance = pow(aiPos.x-otherPos.x, 2) + pow(aiPos.y-otherPos.y, 2) + pow(aiPos.z-otherPos.z, 2);

    if (distance <= m_proximity) return true;
    return false;
}

// at every update, check if AIPos is near otherPos
Status ProximityCondition::update(float seconds){
    // while entity is still pathfinding, keep returning success
    if (m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue) return Status::SUCCESS;
    if (checkProximity() && m_global_blackboard["player"].conditionData["isJumping"].conditionTrue) return Status::SUCCESS;

    return Status::FAIL;
}


void ProximityCondition::reset(){}
void ProximityCondition::addChildren(BTNode *node){}
