#ifndef GameWorld_H
#define GameWorld_H
#include "Game/Environment/environmentsystem.h"
#include "Game/Systems/CollisionSystems/environmentcollisiondetectionsystem.h"
#include "Game/Systems/aisystem.h"
#include "Game/Systems/system.h"
#include "Graphics/global.h"
#include "Graphics/camera.h"
#include <GLFW/glfw3.h>
#include "screen.h"
#include <set>

struct Input {
    // key presses
    int inputVal = 0;

    // mouse clicks
    int xpos;
    int ypos;

    // mouse scrolls
    double distance;

    glm::vec3 translationDir = glm::vec3(0.f);

    bool isActive = false;
    bool isHeld = false;

    bool isClicked = false;
    int checkClickTime = 0;
};


class GameWorld
{
public:
    GameWorld(std::map<int, Input>& input_map);


    void update(double deltaTime);
    void draw();
    void keyEvent(int key, int action);
    void mousePosEvent(double xpos, double ypos);
    void mouseButtonEvent(int button, int action);
    void scrollEvent(double distance);
    void windowResizeEvent(int width, int height);
    void framebufferResizeEvent(int width, int height);


    template <typename T>
    void addSystem(std::shared_ptr<T> &&component){
        m_systems.put<T>(std::forward<std::shared_ptr<T>>(component));
    }


    template <class T>
    T* getSystem(){
        auto comp = m_systems.find<T>();
        assert(comp != m_systems.end());
        return static_cast<T*>(comp->second.get());
    }

    template <class T>
    void removeSystem(){
        m_external_systems.remove<T>();
    }

    template <typename T>
    void addExternalSystem(std::shared_ptr<T> &&component){
        m_external_systems.put<T>(std::forward<std::shared_ptr<T>>(component));
    }


    template <class T>
    T* getExternalSystem(){
        auto comp = m_external_systems.find<T>();
        assert(comp != m_external_systems.end());
        return static_cast<T*>(comp->second.get());
    }

    template <class T>
    void removeExternalSystem(){
        m_external_systems.remove<T>();
    }



private:
//    GLuint VAO, VBO;
//    GLuint grass_tex, wind_tex;
//    std::vector<glm::vec3> positions;

    std::map<int, Input> m_input_map;

    void createKeyInput(int inputVal);
    void initializeInputMap();
    void reset();


    void initializeGrassVBO();

    TypeMap<std::shared_ptr<System>> m_systems;
    TypeMap<std::shared_ptr<System>> m_external_systems;

    std::shared_ptr<Camera> camera;

    std::map<std::string, std::shared_ptr<GameObject>> m_gameobjects;
    std::vector<std::string> m_game_object_names;
    std::map<std::string, std::shared_ptr<GameObject>> m_dynamic_gameobjects;
    std::vector<std::string> m_dynamic_game_object_names;
    std::map<std::string, std::shared_ptr<GameObject>> m_rigid_gameobjects;
    std::vector<std::string> m_rigid_game_object_names;


    std::shared_ptr<EnvironmentCollisionDetectionSystem> collisionSystem;
    std::shared_ptr<EnvironmentSystem> environmentSystem;


    // map that contains submaps of each data type
    std::map<std::string, BlackboardData> m_global_blackboard;
    bool m_checkingMouseClick = false;
    int m_mouseIsHeldTime = 0;

    glm::vec2 m_mousePos;

    std::map<std::string, std::vector<std::shared_ptr<GameObject>>> m_lootables;
    std::set<std::string> m_shownScreens;


};

#endif // GameWorld_H
