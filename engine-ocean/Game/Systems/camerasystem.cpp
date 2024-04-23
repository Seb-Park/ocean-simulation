#include "camerasystem.h"
#include "Game/Components/TransformComponent.h"
#include "Game/GameWorld.h"
#include <memory>
#include <math.h>


CameraSystem::CameraSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                           std::shared_ptr<Camera>& camera,
                           std::map<int, Input>& input_map):
    m_gameobjects(gameobjects),
    m_camera(camera),
    m_input_map(input_map),
    first_person(false),
    third_person(true),
    m_player(m_gameobjects["player"])
{
    camera->setPos(glm::vec3(0.f, 15.f, 0.f));
    camera->setLook(glm::vec3(-28.6761,5.40863,46.5823));


}

glm::vec3 CameraSystem::calculateIdealOffset(glm::vec3 initial){
    //glm::quat playerRotation = glm::quat(
//    float y_rotate = m_player->getComponent<TransformComponent>()->getYRotationAngle();
//    glm::quat playerRot = glm::angleAxis(y_rotate, glm::vec3(0.f, 1.f, 0.f));
    glm::vec4 idealOffset = glm::vec4(initial, 0.f);

    glm::vec4 rotatedOffset = m_player->getComponent<TransformComponent>()->getMT()->getRotation() * idealOffset * .2f;
    glm::vec3 rotated = glm::vec3(rotatedOffset) + m_player->getComponent<TransformComponent>()->getPos();

    return rotated;
}

void CameraSystem::update(double deltaTime){
    if (m_gameobjects["player"]->getComponent<TransformComponent>()->movingLaterally){
        glm::vec3 player_pos = m_gameobjects["player"]->getComponent<TransformComponent>()->getPos();


        glm::vec4 movement = m_distance*glm::vec4(1.f,0.f,0.f, 0.f);
        glm::vec4 transl = m_gameobjects["player"]->getComponent<TransformComponent>()->getMT()->getRotation()*movement;
        transl.y = 0.f;

        glm::vec3 offset = glm::vec3(0.f, 5.f, 0.f);
        glm::vec3 newPos = player_pos + glm::vec3(transl) + offset;

       // arcballRotation(newPos.x, newPos.z);
        m_camera->setPos(newPos);

        if (glm::cross(m_camera->getUp(), player_pos) != glm::vec3(0.f)){
            m_camera->setLook(player_pos);
        }
    }
}


void CameraSystem::draw(){
    Global::graphics.bindShader("phong");
   Global::graphics.setCameraData(m_camera);
}

void CameraSystem::firstPersonMode(){
    first_person = true;
    third_person = false;
    m_camera->setPos(m_gameobjects["player"]->getComponent<TransformComponent>()->getPos());
}

void CameraSystem::thirdPersonMode(glm::vec3 curr_cam_pos){
    first_person = false;
    third_person = true;
    //camera->setPos(curr_cam_pos);
}

void CameraSystem::scrollEvent(double distance){
   // calculateZoom(distance);
//    m_distance += distance;
//    if (m_distance < 5.f) m_distance = 5.f;
//    if (m_distance > 40.f) m_distance = 40.f;
}

void CameraSystem::mousePosEvent(double xpos, double ypos){

    // if right mouse button is pressed, calculate change of the new input x and y pos
    if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].isHeld){
        float deltaX = xpos - prev_mouse_pos.x;
        float deltaY = ypos - prev_mouse_pos.y;

        // upate prev mouse pos
        prev_mouse_pos = glm::vec2(xpos, ypos);

        arcballRotation(deltaX, deltaY);

    }
}

void CameraSystem::calculateZoom(float scrollDistance){
    float zoomLevel = scrollDistance*.1f;
    m_distanceFromPlayer -= zoomLevel;
}

void CameraSystem::calculatePitch(double deltaY){
    if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].isHeld){
        float pitchChange = deltaY*.1f;
        m_pitch -= pitchChange;
    }
}

void CameraSystem::calculateAngleAroundPlayer(double deltaX){
    if (m_input_map[GLFW_MOUSE_BUTTON_LEFT].isHeld){
        float angleChange = deltaX*.3f;
        m_angleAroundPlayer -= angleChange;
    }
}

void CameraSystem::arcballRotation(float deltaX, float deltaY){
    glm::vec4 cam_pos = glm::vec4(m_camera->getPos(), 1.f);
    glm::vec4 pivot = glm::vec4(m_player->getComponent<TransformComponent>()->getPos(), 1.f);

    // calculate rotation amount given mouse movement
    float deltaAngleX = (2*M_PI / static_cast<float>(m_camera->getWidth()));
    float deltaAngleY = (M_PI / static_cast<float>(m_camera->getHeight()));

    float xAngle = deltaX * deltaAngleX;
    float yAngle = deltaY * deltaAngleY;

    // handles when camera dir is same as up vector
    float cosAngle = glm::dot(m_camera->getViewDirection(), m_camera->getUp());
    float sign = 1.f;
    if (yAngle < 0) sign = -1.f;
    if (cosAngle * sign > .99f) yAngle = 0;

    // rotate camera around pivtor point
    glm::mat4 rotX(1.f);
    rotX = glm::rotate(rotX, xAngle, m_camera->getUp());
    cam_pos = (rotX * (cam_pos - pivot)) + pivot;

    glm::mat4 rotY(1.f);
    rotY = glm::rotate(rotY, yAngle, m_camera->getRight());
    glm::vec3 final_pos = (rotY * (cam_pos - pivot)) + pivot;

    // update camera view
    m_camera->setPos(final_pos);

}

void CameraSystem::setCameraPos(){
//    float horiz_distance = m_distanceFromPlayer*glm::cos(glm::radians(m_pitch));
//    float vert_distance = m_distanceFromPlayer*glm::sin(glm::radians(m_pitch));

//    glm::vec3 player_pos = m_gameobjects["player"]->getComponent<TransformComponent>()->getPos();
//    float y_camPos = player_pos.y + vert_distance;

//    float theta = glm::radians(m_angleAroundPlayer);
//    float offsetX = horiz_distance*glm::cos((theta));
//    float offsetZ = horiz_distance*glm::sin((theta));

//    glm::vec3 camPos = glm::vec3(player_pos.x - offsetX, y_camPos, player_pos.z - offsetZ);

   // m_yaw = glm::radians(180-theta);

//    glm::vec3 lookVector = camPos +
//            glm::vec3(glm::cos(pitch)*glm::sin(m_yaw), glm::sin(pitch), glm::cos(pitch)*glm::cos(m_yaw));


//    glm::vec4 distanceCamToPlayer = glm::vec4(camPos - player_pos, 0.f);
//    glm::vec4 look =  m_gameobjects["player"]->getComponent<TransformComponent>()->getMT()->getRotation() * distanceCamToPlayer;

//    //glm::vec3 idealOffset = (glm::vec3(-7.f, 5.f, -10.f)) + m_player->getComponent<TransformComponent>()->getPos();
//    glm::vec3 idealLookat = calculateIdealOffset(glm::vec3(0.f, -2.f, 40.f));

//    m_currPos = idealOffset;
//    m_currLook = idealLookat;

//    glm::vec3 direction;
//    direction.x = cos(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
//    direction.y = sin(glm::radians(m_pitch));
//    direction.z = sin(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));

//    // calculate look
//    glm::vec4 rotatedPlayerPoint = m_player->getComponent<TransformComponent>()->getMT()->getRotation()* glm::vec4(player_pos, 1.f);
//    float distanceX = glm::distance(rotatedPlayerPoint.x, player_pos.x);
//    float distanceY = glm::distance(rotatedPlayerPoint.y, player_pos.y);
//    float distanceZ = glm::distance(rotatedPlayerPoint.z, player_pos.z);



    //camPos = camPos + glm::vec3(distanceX, distanceY, distanceZ);
    //glm::vec3 newLook = camPos - player_pos;



   // glm::vec4 lookToPlayer = glm::vec4(camPos - player_pos, 0.f);
   // glm::vec4 rotatedLook = m_player->getComponent<TransformComponent>()->getMT()->getRotation() * lookToPlayer;
   // float distLookToPlayer = glm::distance(glm::vec3(rotatedLook), player_pos);


   // m_camera->setPos(camPos);
    //m_camera->translate(m_currPos);
    //m_camera->rotate(m_yaw, glm::vec3(0.f, 1.f, 0.f));
    //m_camera->setLook(newLook);


   // m_camera->setPos(glm::vec3(x_camPos, y_camPos, z_camPos));




}

