#ifndef PHYSICSSYSTEM_H
#define PHYSICSSYSTEM_H
#include "Game/Components/TransformComponent.h"
#include "system.h"


class PhysicsSystem : public System
{
public:
    PhysicsSystem(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                  std::map<std::string, BlackboardData>& global_blackboard,
                  std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;

private:
    float gravitySimulation(float &initial_v, double deltaTime, float snapshot_time, float gravity);
    TransformComponent* getTransform(std::shared_ptr<GameObject> &go);

    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::map<std::string, BlackboardData>& m_global_blackboard;
    std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& m_lootables;

    float horiz_velocity = .005f;
    float snapshot_time = 0.f;
};

#endif // PHYSICSSYSTEM_H
