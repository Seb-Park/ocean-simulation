#ifndef ENVIRONMENTCOLLISIONDETECTIONSYSTEM_H
#define ENVIRONMENTCOLLISIONDETECTIONSYSTEM_H


#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include <memory>
#include "Game/Systems/system.h"

class EnvironmentCollisionDetectionSystem /*: public System*/
{
public:
    EnvironmentCollisionDetectionSystem(
            std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
            std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
            std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables,
            std::map<std::string, BlackboardData>& global_blackboard);
    void draw() ;
    void update(double deltaTime);
    void scrollEvent(double distance);
    void mousePosEvent(double xpos, double ypos) ;


private:
    void detectCollisionWithEnvironment(double deltaTime);
    CollisionComponent* getCollisionComp(std::shared_ptr<GameObject> &go);
    TransformComponent* getTransform(std::shared_ptr<GameObject> &go);

    void collideDynamic(std::string go_name, std::shared_ptr<GameObject> go);


    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::unique_ptr<EllipsoidTriangleCollisionSystem> m_ellipsoid_triangle_collision_system;
    //std::unique_ptr<BVHTree> m_bvh_system;
    std::vector<std::shared_ptr<GameObject>> all_rigid_go;
    std::map<std::string, BlackboardData>& m_global_blackboard;

    std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& m_lootables;

};

#endif // ENVIRONMENTCOLLISIONDETECTIONSYSTEM_H
