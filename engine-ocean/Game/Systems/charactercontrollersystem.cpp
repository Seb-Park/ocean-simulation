#include "charactercontrollersystem.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include "Game/Components/TransformComponent.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include <glm/gtx/quaternion.hpp>

CharacterControllerSystem::CharacterControllerSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                                                     std::shared_ptr<Camera>& camera_param,
                                                     std::map<int, Input>& input_map,
                                                     std::map<std::string, BlackboardData>& global_blackboard):
    m_gameobjects(gameobjects),
    camera(camera_param),
    m_input_map(input_map),
    m_global_blackboard(global_blackboard)
{
}

void CharacterControllerSystem::draw(){
}

TransformComponent* CharacterControllerSystem::getPlayerTransform(){
    return m_gameobjects.at("player")->getComponent<TransformComponent>();
}

glm::vec3 CharacterControllerSystem::getPlayerPos(){
    return m_gameobjects.at("player")->getComponent<TransformComponent>()->getPos();
}


float CharacterControllerSystem::jumpPlayer(float &initial_v, double deltaTime, float snapshot_time, float gravity){
    float t = deltaTime-snapshot_time;
    float delta_y = initial_v*t + (.5f)*gravity*t*t;
    return delta_y;
}

bool CharacterControllerSystem::movePlayerLaterally(glm::vec3 dir, glm::vec3 perp, glm::vec3 &m_pos, float dt){
    glm::vec3 translationDir;
    bool moving_laterally = false;
    std::shared_ptr<ModelTransform> temp_mt = getPlayerTransform()->getMT();

    if (m_input_map[GLFW_KEY_W].isActive){
        m_currentSpeed = RUN_SPEED;
        moving_laterally = true;
    } else if (m_input_map[GLFW_KEY_S].isActive){
        m_currentSpeed = -RUN_SPEED;
        moving_laterally = true;
    } else {
        m_currentSpeed = 0.f;
    }


    glm::mat4 rot = glm::mat4(1.f);
    // player turning
    if (m_input_map[GLFW_KEY_D].isActive){
        m_currentTurnSpeed = -TURN_SPEED;
        moving_laterally = true;
    } else if (m_input_map[GLFW_KEY_A].isActive){
        m_currentTurnSpeed = TURN_SPEED;
        moving_laterally = true;
    } else {
        m_currentTurnSpeed = 0.f;
    }

    // rotate player
    getPlayerTransform()->getMT()->rotate(m_currentTurnSpeed*dt, glm::vec3(0.f, 1.f, 0.f));


    //translate player based on rotation
    float distance = m_currentSpeed*dt;
    glm::vec4 movement = distance*glm::vec4(1.f,0.f,0.f, 0.f);
    glm::vec4 transl = getPlayerTransform()->getMT()->getRotation()*movement;

    translationDir = glm::vec3(transl.x, 0.f, transl.z);
    temp_mt->translate(translationDir);
    glm::vec3 potential_pos = temp_mt->getPos();
    m_pos = potential_pos;


    return moving_laterally;
}

void CharacterControllerSystem::handlePlayerMovement(double deltaTime){
    // timing
        float dt = deltaTime - snapshot_time;
        snapshot_time = deltaTime;

        // get player positions
        m_pos = getPlayerPos();
        getPlayerTransform()->old_pos = m_pos;

        // movement
        getPlayerTransform()->movingLaterally = false;
        float speedFactor = horiz_velocity * deltaTime;
        glm::vec3 look = camera->getLook();
        glm::vec3 dir = glm::normalize(glm::vec3(look.x, 0.f, look.z)) * speedFactor;
        glm::vec3 perp = glm::vec3(look.z, 0.f, -look.x) * speedFactor;

        // updates m_pos
        if (movePlayerLaterally(dir, perp, m_pos, dt)){
            getPlayerTransform()->movingLaterally = true;
        }

        if (m_input_map[GLFW_KEY_SPACE].isActive && (getPlayerTransform()->onGround)){
            getPlayerTransform()->gravity = -25.f;
            getPlayerTransform()->yVelocity = 20.f;
            getPlayerTransform()->onGround = false;
            //if (!gravity) gravity = true;
            m_global_blackboard["player"].conditionData["isJumping"].conditionTrue = true;
        } else if (!m_input_map[GLFW_KEY_SPACE].isActive && getPlayerTransform()->onGround){
            m_global_blackboard["player"].conditionData["isJumping"].conditionTrue = false;
        }

        // PLAYER-SPECIFIC GRAVITY (comment these two lines to keep player still)
        //if (gravity){
        m_pos.y += jumpPlayer(getPlayerTransform()->yVelocity, dt, 0, getPlayerTransform()->gravity);
        getPlayerTransform()->yVelocity = getPlayerTransform()->yVelocity + getPlayerTransform()->gravity*(dt);
        //}

        // store m_pos as estimated final_pos
        m_global_blackboard["player"].locationData.setToPos = m_pos;

        getPlayerTransform()->estimated_final_pos = m_pos;
}


void CharacterControllerSystem::update(double deltaTime){
    handlePlayerMovement(deltaTime);
}

void CharacterControllerSystem::scrollEvent(double distance){}
void CharacterControllerSystem::mousePosEvent(double xpos, double ypos){}
