#include "environmentcollisiondetectionsystem.h"
#include "Game/Components/CollisionComponents/BoundingDynamicMesh.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include <utility>

EnvironmentCollisionDetectionSystem::EnvironmentCollisionDetectionSystem(
        std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
        std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
        std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables,
        std::map<std::string, BlackboardData>& global_blackboard) :
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_ellipsoid_triangle_collision_system(std::make_unique<EllipsoidTriangleCollisionSystem>(rigid_gameobjects)),
    m_global_blackboard(global_blackboard),
    m_lootables(lootables)
{

}

TransformComponent* EnvironmentCollisionDetectionSystem::getTransform(std::shared_ptr<GameObject> &go){
    return go->getComponent<TransformComponent>();

}

CollisionComponent* EnvironmentCollisionDetectionSystem::getCollisionComp(std::shared_ptr<GameObject> &go){
    return go->getComponent<CollisionComponent>();
}

void EnvironmentCollisionDetectionSystem::detectCollisionWithEnvironment(double deltaTime){
    for (auto &go : m_dynamic_gameobjects){
        collideDynamic(go.first, go.second);

    }
}


void EnvironmentCollisionDetectionSystem::collideDynamic(std::string go_name, std::shared_ptr<GameObject> go){
    if (!go->hasComponent<CollisionComponent>()){
        return;
    }

    std::pair<std::vector<CollisionData>, glm::vec3> pair =
        m_ellipsoid_triangle_collision_system->mtvSlide(getTransform(go)->old_pos,
                                                        //getTransform(go.second)->estimated_final_pos,
                                                        m_global_blackboard[go_name].locationData.setToPos,
                                                        getCollisionComp(go)->getCollisionShape<BoundingDynamicMesh>()->getEllipsoidDimensions());

    // assume not on ground first; then determine if we are touching ground
    getTransform(go)->onGround = false;
    m_global_blackboard[go_name].conditionData["onGround"].conditionTrue = false;
    for (const CollisionData &collision : pair.first){
        if (glm::dot(glm::vec3(0,1,0), collision.triangle_n) > 0.f){
           getTransform(go)->onGround = true;
           m_global_blackboard[go_name].conditionData["onGround"].conditionTrue = true;

        }
    }

    if (getTransform(go)->onGround){
        if (getTransform(go)->yVelocity < 0){
            getTransform(go)->yVelocity = 0.f;
            getTransform(go)->gravity = 0.f;
        }
    }

    if (!getTransform(go)->onGround){
        getTransform(go)->gravity = -25.f;
    }


    if (!getTransform(go)->movingLaterally &&
            getTransform(go)->onGround){
        getTransform(go)->setPos(getTransform(go)->old_pos);
    } else {
        getTransform(go)->setPos(pair.second);
    }

    getCollisionComp(go)->getCollisionShape<BoundingDynamicMesh>()->updateCenterPos(getTransform(go)->getPos());
}


void EnvironmentCollisionDetectionSystem::update(double deltaTime){
    detectCollisionWithEnvironment(deltaTime);
}

void EnvironmentCollisionDetectionSystem::draw(){}
void EnvironmentCollisionDetectionSystem::scrollEvent(double distance){}
void EnvironmentCollisionDetectionSystem::mousePosEvent(double xpos, double ypos){}



