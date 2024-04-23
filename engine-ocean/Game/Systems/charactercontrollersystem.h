#ifndef CHARACTERCONTROLLERSYSTEM_H
#define CHARACTERCONTROLLERSYSTEM_H


#include "Game/Components/TransformComponent.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/GameWorld.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include <map>



class CharacterControllerSystem : public System
{
public:
    CharacterControllerSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                              std::shared_ptr<Camera>& camera_param,
                              std::map<int, Input>& input_map,
                              std::map<std::string, BlackboardData>& global_blackboard);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;

private:
    bool movePlayerLaterally(glm::vec3 dir, glm::vec3 perp, glm::vec3 &m_pos, float dt);
    void handlePlayerMovement(double deltaTime);
    TransformComponent* getPlayerTransform();

    std::map<std::string, std::shared_ptr<GameObject>>& m_gameobjects;
    std::shared_ptr<Camera>& camera;
    std::map<int, Input>& m_input_map;

    void movePlayer(glm::vec3 direction, int key);
    void setPlayerPos(glm::vec3 new_pos);
    glm::vec3 getPlayerPos();
    float jumpPlayer(float &initial_v, double deltaTime, float snapshot_time, float gravity);


    glm::vec3 m_pos = glm::vec3(0.f);
    float horiz_velocity = .005f;
    float snapshot_time = 0.f;

    std::map<std::string, BlackboardData>& m_global_blackboard;

    float m_currentSpeed = 0.f;
    float m_currentTurnSpeed = 0.f;

    float RUN_SPEED = 8.f;
    float TURN_SPEED = 2.f;

    glm::quat m_quat;
    glm::vec3 m_axis;

    float HALFPI = 0.78539f;
    bool gravity = false;



};

#endif // CHARACTERCONTROLLERSYSTEM_H
