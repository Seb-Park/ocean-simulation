#ifndef CAMERASYSTEM_H
#define CAMERASYSTEM_H


#include "Game/GameWorld.h"
#include "Graphics/camera.h"
#include <memory>
#include "system.h"

class CameraSystem : public System
{
public:
    CameraSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                 std::shared_ptr<Camera>& camera,
                 std::map<int, Input>& input_map);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;

private:
    void firstPersonMode();
    void thirdPersonMode(glm::vec3 curr_cam_pos);

    void setCameraPos();
    void calculateAngleAroundPlayer(double deltaX);
    void calculatePitch(double deltaY);
    void calculateZoom(float scrollDistance);

    glm::vec3 calculateIdealOffset(glm::vec3 initial);
    void arcballRotation(float deltaX, float deltaY);
    float m_distance = -10.f;










    std::shared_ptr<Camera>& m_camera;
    std::map<std::string, std::shared_ptr<GameObject>>& m_gameobjects;
    std::map<int, Input>& m_input_map;
    bool first_person = false;
    bool third_person = true;
    float first_person_offset = 0.01f;
    glm::vec2 prev_mouse_pos = glm::vec2(0.f);
    float horiz_velocity = .005f;

    std::shared_ptr<GameObject> m_player;
    float m_distanceFromPlayer = 50.f;
    float m_angleAroundPlayer = 0.f;
    float m_pitch = 20.f;
    float m_yaw = 0.f;

    bool m_isInverted = false;
    float m_angle = 0.f;
    float m_height = 0.f;
    int m_invertY = 1;

    glm::vec3 m_currPos;
    glm::vec3 m_currLook;

    float m_currentTurnSpeed = 0.f;

    float TURN_SPEED = 2.f;
    float snapshot_time = 0.f;


};

#endif // CAMERASYSTEM_H
