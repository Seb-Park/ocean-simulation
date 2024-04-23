#include "physicssystem.h"
#include "Game/Components/TransformComponent.h"

PhysicsSystem::PhysicsSystem(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                             std::map<std::string, BlackboardData>& global_blackboard,
                             std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables) :
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_global_blackboard(global_blackboard),
    m_lootables(lootables)

{

}

TransformComponent* PhysicsSystem::getTransform(std::shared_ptr<GameObject> &go){
    return go->getComponent<TransformComponent>();
}

float PhysicsSystem::gravitySimulation(float &initial_v, double deltaTime, float snapshot_time, float gravity){
    float t = deltaTime-snapshot_time;
    float delta_y = initial_v*t + (.5f)*gravity*t*t;
    return delta_y;
}


void PhysicsSystem::update(double deltaTime){
    //std::cout << "physics" << std::endl;
    float dt = deltaTime - snapshot_time;
    snapshot_time = deltaTime;

    for (auto &go : m_dynamic_gameobjects){
        // player has its own physics
        if (go.first != "player"){
            // position
            glm::vec3 m_pos = getTransform(go.second)->getPos();
            getTransform(go.second)->old_pos = m_pos;

            // effect it by gravity
            m_pos.y += gravitySimulation(getTransform(go.second)->yVelocity, dt, 0, getTransform(go.second)->gravity);
            getTransform(go.second)->yVelocity = getTransform(go.second)->yVelocity + getTransform(go.second)->gravity*(dt);

            m_global_blackboard[go.first].locationData.setToPos = m_pos;
            // store m_pos as estimated final_pos
            getTransform(go.second)->estimated_final_pos = m_pos;
        }
    }

//    for (auto &lootGroup : m_lootables){
//        for (auto &loot : lootGroup.second){
//            // position
//            glm::vec3 m_pos = getTransform(loot)->getPos();
//            getTransform(loot)->old_pos = m_pos;

//            // effect it by gravity
//            m_pos.y += gravitySimulation(getTransform(loot)->yVelocity, dt, 0, getTransform(loot)->gravity);
//            getTransform(loot)->yVelocity = getTransform(loot)->yVelocity + getTransform(loot)->gravity*(dt);

//            // store m_pos as estimated final_pos
//            getTransform(loot)->estimated_final_pos = m_pos;
//        }

//    }

}




void PhysicsSystem::draw(){}
void PhysicsSystem::scrollEvent(double distance){}
void PhysicsSystem::mousePosEvent(double xpos, double ypos){}

