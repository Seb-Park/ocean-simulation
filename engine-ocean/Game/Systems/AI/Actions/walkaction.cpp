#include "walkaction.h"
#include "Game/Components/PathfindComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include "glm/glm.hpp"
#include <memory>


WalkAction::WalkAction(std::string entity_id,
                       std::map<std::string, BlackboardData>& global_blackboard):
    m_global_blackboard(global_blackboard)
{

    m_entity_id = entity_id;
    m_path.clear();

    m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue = false;
    m_global_blackboard[m_entity_id].conditionData["atDestination"].conditionTrue = false;
    m_global_blackboard[m_entity_id].conditionData["pathfound"].conditionTrue = false;




}

void WalkAction::setPath(glm::vec3 entity_pos){
    std::cout << "---------SETTING PATHH" << std::endl;
    m_destination = m_global_blackboard["player"].locationData.currPos;
    m_path = m_global_blackboard["navmesh"].environment->getComponent<PathfindComponent>()->getPath(glm::vec3(-0.58249, 0, -0.0210782), glm::vec3(19.5371, 0, 1.39167));
}

// only activates if the previous conditions are true
Status WalkAction::update(float seconds){

    glm::vec3 pos =  m_global_blackboard[m_entity_id].locationData.currPos;

    // get a path if entity is not pathfinding
    if (!m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue &&
            !m_global_blackboard[m_entity_id].conditionData["atDestination"].conditionTrue &&
            !m_global_blackboard[m_entity_id].conditionData["pathfound"].conditionTrue){
        setPath(pos);
        m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue = true;
        m_global_blackboard[m_entity_id].conditionData["pathfound"].conditionTrue = true;

    }

    if (!m_path.empty()){
        if (m_global_blackboard[m_entity_id].conditionData["onGround"].conditionTrue){
            //std::cout << "on ground" << std::endl;

            std::shared_ptr<ModelTransform> temp_mt = std::make_shared<ModelTransform>();
            temp_mt->setPos(pos);
            glm::vec3 v = m_path.back();//glm::vec3(m_path.back().x, pos.y, m_path.back().z);
            glm::vec3 dir = glm::normalize(v-temp_mt->getPos());
            temp_mt->translate(dir);
            //std::cout << "v: (" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
            glm::vec3 pos_eps = v + .01f;
            glm::vec3 neg_eps = v - .01f;

            // pop if entity within a certain episilon of node
            if (neg_eps.x < temp_mt->getPos().x < pos_eps.x &&
                    neg_eps.z < temp_mt->getPos().z < pos_eps.z){
                 m_path.pop_back();
            }

            m_global_blackboard[m_entity_id].locationData.setToPos = temp_mt->getPos();
        }

        return Status::RUNNING;
    }

    // if reached destination, then walking succeeded
    if (pos.x == m_destination.x && pos.z == m_destination.z){
        m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue = false;
        m_global_blackboard[m_entity_id].conditionData["atDestination"].conditionTrue = true;
        std::cout << "-reached-" << std::endl;
        return Status::SUCCESS;

    }

    // otherwise
    m_global_blackboard[m_entity_id].conditionData["isPathfinding"].conditionTrue = false;
    return Status::FAIL;

}

void WalkAction::reset(){}
void WalkAction::addChildren(BTNode *node){}


