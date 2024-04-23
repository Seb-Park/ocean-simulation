#ifndef COLLISIONSYSTEM_H
#define COLLISIONSYSTEM_H
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameWorld.h"
#include "Game/Systems/CollisionSystems/UniformGrid/uniformgrid.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include "Game/Systems/system.h"

class CollisionSystem : public System
{
public:
    CollisionSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                                     std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                     std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                     std::map<int, Input>& input_map,
                                     std::map<std::string, BlackboardData>& global_blackboard);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;
private:
    bool areCollidable(const std::shared_ptr<GameObject> &a, const std::shared_ptr<GameObject> &b);
    void detectCylinderCollisions(std::shared_ptr<GameObject> &a, std::shared_ptr<GameObject> &b, float deltaTime,
                                  std::string a_name, std::string b_name);
    void resolveCollision(std::shared_ptr<GameObject> &go, glm::vec3 translation, float deltaTime, std::string go_name);
    bool checkCollisionShape(const std::shared_ptr<GameObject> &go, std::string shape_type);
    glm::vec3 projectMotion(const glm::vec3 &initial_pos, float deltaTime);

    TransformComponent* getTransform(std::shared_ptr<GameObject> &go);
    CollisionComponent* getCollisionComp(std::shared_ptr<GameObject> &go);

    void rigidDynamicCollisions(double deltaTime);
    void dynamicDynamicCollisions(double deltaTime);


    // cylinder-cylinder
    glm::vec3 collideCylinderCylinder(const std::shared_ptr<GameObject> &a_go, const std::shared_ptr<GameObject> &b_go);
    float calculateLineMVT(const Cylinder &a, const Cylinder &b);
    glm::vec2 calculateCircleMVT(const Cylinder &a, const Cylinder &b);


    std::map<std::string, std::shared_ptr<GameObject>>& m_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;

    std::unique_ptr<UniformGrid> m_uniform_grid_system;
    std::map<std::string, BlackboardData>& m_global_blackboard;
};


#endif // COLLISIONSYSTEM_H
