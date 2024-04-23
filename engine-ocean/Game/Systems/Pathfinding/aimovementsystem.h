#ifndef AIMOVEMENTSYSTEM_H
#define AIMOVEMENTSYSTEM_H
#include "Game/Components/TransformComponent.h"
#include "Game/Systems/Pathfinding/pathfinder.h"
#include "Game/Systems/system.h"


class AIMovementSystem : public System
{
public:
    AIMovementSystem(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                       std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;

private:
    float gravitySimulation(float &initial_v, double deltaTime, float snapshot_time, float gravity);
    TransformComponent* getTransform(std::shared_ptr<GameObject> &go);

    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;


    float horiz_velocity = .005f;
    float snapshot_time = 0.f;

    std::vector<glm::vec3> m_path;

    //std::shared_ptr<Pathfinder> m_pathfinder;
};

#endif // AIMOVEMENTSYSTEM_H
