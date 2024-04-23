#ifndef OBJECTCREATIONSYSTEM_H
#define OBJECTCREATIONSYSTEM_H


#include "Game/GameObjects/GameObject.h"
#include "Game/Systems/aisystem.h"
#include "Game/Ocean/ocean.h"

#include "system.h"


class ObjectCreationSystem : public System
{
public:
    ObjectCreationSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                                               std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                               std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                               std::map<std::string, BlackboardData>& global_blackboard,
                         std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables);
    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;

private:
    void initializeAllObjects();
    void initializePlayerObject();
    void initializeGroundObject();
    void initializeObstacle();
    void generateBamboo();
    void generateBall(glm::vec3 scale, glm::vec3 pos, int number);
    void addObject_NoTransform(std::string shape_name, std::string obj_name, std::string filename);
    void makeBVHDemo();
    void initializeBackground();
    void initOcean();



    void makeLootable(std::string shapename, std::string filename, int count, glm::vec3 scale);





    void initializeCeilingMesh();
    void initializeSlopedGround();
    void initializeGround();
    void addLight();
    void makeNavMesh();




    void insertAnyObject(const std::string name, const std::shared_ptr<GameObject> &game_obj);
    void insertDynamicObject(const std::string name, const std::shared_ptr<GameObject> &game_obj);
    void insertRigidObject(const std::string name, const std::shared_ptr<GameObject> &game_obj);




    std::map<std::string, std::shared_ptr<GameObject>>& m_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;

    std::map<std::string, BlackboardData>& m_global_blackboard;

    std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& m_lootables;


    float m_groundBounds = 0.f;
    float m_ground_level = 0.f;


    ocean m_ocean;
    std::shared_ptr<GameObject> m_ocean_shape;
    double m_time = 0.0;
    double m_timestep = .01;


};

#endif // OBJECTCREATIONSYSTEM_H
