#include "aimovementsystem.h"
#include "Game/Components/TransformComponent.h"
#include "Game/Components/PathfindComponent.h"
#include "Game/Systems/AI/aibehaviorcomponent.h"

AIMovementSystem::AIMovementSystem(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                   std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects) :
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_rigid_gameobjects(rigid_gameobjects)

{
   // m_path = m_rigid_gameobjects.at("navmesh")->getComponent<PathfindComponent>()->getPath(glm::vec3(-0.58249, 0, -0.0210782), glm::vec3(19.5371, 0, 1.39167));

}

TransformComponent* AIMovementSystem::getTransform(std::shared_ptr<GameObject> &go){
    return go->getComponent<TransformComponent>();
}


void AIMovementSystem::update(double deltaTime){
    for (auto &go : m_dynamic_gameobjects){
        if (go.second->hasComponent<AIBehaviorComponent>()){
            go.second->getComponent<AIBehaviorComponent>()->update(deltaTime);
        }
    }
}




void AIMovementSystem::draw(){}
void AIMovementSystem::scrollEvent(double distance){}
void AIMovementSystem::mousePosEvent(double xpos, double ypos){}

