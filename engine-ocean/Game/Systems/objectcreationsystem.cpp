#include "objectcreationsystem.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/boundingtriangle.h"
#include "Game/Components/DrawComponent.h"
#include "Game/Components/PathfindComponent.h"
#include "Game/Components/TransformComponent.h"
#include "Game/Ocean/ocean.h"

#include <random>

ObjectCreationSystem::ObjectCreationSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                                           std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                           std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                           std::map<std::string, BlackboardData>& global_blackboard,
                                           std::map<std::string, std::vector<std::shared_ptr<GameObject>>>& lootables):
    m_gameobjects(gameobjects),
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_rigid_gameobjects(rigid_gameobjects),
    m_global_blackboard(global_blackboard),
    m_lootables(lootables)

{
    initializeAllObjects();
}

void ObjectCreationSystem::initializeAllObjects(){

    m_ground_level = -.5f;
    initializePlayerObject();
   // initializeSlopedGround();
//    initializeGround();
//    initializeBackground();
    initOcean();

    addLight();
}

void ObjectCreationSystem::initializeSlopedGround(){
    std::shared_ptr<GameObject> sloped_ground = std::make_shared<GameObject>();
    std::vector<glm::vec3> obj_data = Global::graphics.addShape("sloped_ground", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/floor.obj");

    std::shared_ptr<ModelTransform> mt = std::make_shared<ModelTransform>();
    mt->setScale(1.f);
    mt->setPos(glm::vec3(0.f));

    sloped_ground->addComponent<DrawComponent>(std::make_unique<DrawComponent>(Global::graphics.getShape("sloped_ground")));
    sloped_ground->getComponent<DrawComponent>()->addMaterial("grass_tedxxx", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/mossyground.png");
    sloped_ground->addComponent<TransformComponent>(std::make_unique<TransformComponent>(mt, "sloped_ground", m_global_blackboard));
    sloped_ground->addComponent<CollisionComponent>(std::make_unique<CollisionComponent>("obj", obj_data, mt));

    insertRigidObject("sloped_ground", sloped_ground);
}

void ObjectCreationSystem::insertAnyObject(const std::string name, const std::shared_ptr<GameObject> &game_obj){
    m_gameobjects.insert(std::pair<const std::string, std::shared_ptr<GameObject>>(name, game_obj));
}

void ObjectCreationSystem::insertRigidObject(const std::string name, const std::shared_ptr<GameObject> &game_obj){
    m_rigid_gameobjects.insert(std::pair<const std::string, std::shared_ptr<GameObject>>(name, game_obj));
    insertAnyObject(name, game_obj);
}

void ObjectCreationSystem::insertDynamicObject(const std::string name, const std::shared_ptr<GameObject> &game_obj){
    m_dynamic_gameobjects.insert(std::pair<const std::string, std::shared_ptr<GameObject>>(name, game_obj));
    insertAnyObject(name, game_obj);
}

void ObjectCreationSystem::initializePlayerObject(){
    //std::shared_ptr<Shape> shape = Global::graphics.getShape("sphere");
    std::vector<glm::vec3> obj_data = Global::graphics.addShape_withMaterial("mouse", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/mouse2-4.obj",
                                                                             "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/mouse2-4.mtl", true);
    std::shared_ptr<ModelTransform> mt = std::make_shared<ModelTransform>();
    mt->setScale(glm::vec3(.5f));
    mt->setPos(glm::vec3(0.f));

    std::shared_ptr<GameObject> player = std::make_shared<GameObject>();

    player->addComponent<DrawComponent>(std::make_unique<DrawComponent>(Global::graphics.getShapeGroup("mouse")));
    player->addComponent<TransformComponent>(std::make_unique<TransformComponent>(mt, "player", m_global_blackboard, true));
    player->addComponent<CollisionComponent>(std::make_unique<CollisionComponent>("dynamic_mesh", mt, mt->getPos(), obj_data));

    insertDynamicObject("player", player);
}

void ObjectCreationSystem::initializeGround(){
   std::shared_ptr<GameObject> ground = std::make_shared<GameObject>();
   std::vector<glm::vec3> obj_data = Global::graphics.addShape_withMaterial("ground", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/meadow_ground.obj",
                                                              "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/meadow_ground.mtl", true);

    //std::vector<glm::vec3> obj_data = Global::graphics.addShape("ground", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/testplane.obj");
    std::shared_ptr<ModelTransform> mt = std::make_shared<ModelTransform>();
    mt->setPos(glm::vec3(0.f, 0.f, 0.f));
    ground->addComponent<DrawComponent>(std::make_unique<DrawComponent>(Global::graphics.getShapeGroup("ground")));
    ground->addComponent<TransformComponent>(std::make_unique<TransformComponent>(mt, "ground", m_global_blackboard));
    ground->addComponent<CollisionComponent>(std::make_unique<CollisionComponent>("obj", obj_data, mt, true));

    insertRigidObject("ground", ground);
}

void ObjectCreationSystem::initializeBackground(){
    std::shared_ptr<GameObject> bg = std::make_shared<GameObject>();

    // "Snowy Mountain - Terrain" (https://skfb.ly/6RzJV) by artfromheath is licensed under Creative Commons Attribution (http://creativecommons.org/licenses/by/4.0/).
    std::vector<glm::vec3> obj_data = Global::graphics.addShape_withMaterial("bg", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/howl_field_background.obj",
                                                               "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/howl_field_background.mtl", true);

     //std::vector<glm::vec3> obj_data = Global::graphics.addShape("ground", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Meshes/testplane.obj");
     std::shared_ptr<ModelTransform> mt = std::make_shared<ModelTransform>();
     mt->setPos(glm::vec3(0.f, 0.f, 0.f));
     bg->addComponent<DrawComponent>(std::make_unique<DrawComponent>(Global::graphics.getShapeGroup("bg")));
     bg->addComponent<TransformComponent>(std::make_unique<TransformComponent>(mt, "bg", m_global_blackboard));
     bg->addComponent<CollisionComponent>(std::make_unique<CollisionComponent>("obj", obj_data, mt, true));

     insertRigidObject("bg", bg);

}

void ObjectCreationSystem::initOcean(){
    m_ocean_shape = std::make_shared<GameObject>();

    std::vector<glm::vec3> obj_data = Global::graphics.addShape_manual("ocean", m_ocean.get_vertices(), m_ocean.get_faces(), true);


    std::shared_ptr<ModelTransform> mt = std::make_shared<ModelTransform>();
    mt->setScale(1.f);
    mt->setPos(glm::vec3(0.f, 0.f, 0.f));
    m_ocean_shape->addComponent<DrawComponent>(std::make_unique<DrawComponent>(Global::graphics.getShape("ocean")));
    m_ocean_shape->getComponent<DrawComponent>()->addMaterial("grass_tedxxx", "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/mossyground.png");

    m_ocean_shape->addComponent<TransformComponent>(std::make_unique<TransformComponent>(mt, "ocean", m_global_blackboard));
    m_ocean_shape->addComponent<CollisionComponent>(std::make_unique<CollisionComponent>("obj", obj_data, mt, true));

    insertRigidObject("ocean", m_ocean_shape);


}


void ObjectCreationSystem::addLight(){
    std::shared_ptr<Light> light1 = std::make_shared<Light>(LightType::DIRECTIONAL, glm::vec3(0,20,6.5));
    std::vector<std::shared_ptr<Light>> lights;
    lights.push_back(light1);

    Global::graphics.bindShader("phong");
    Global::graphics.setLights(lights);
}


void ObjectCreationSystem::draw(){}
void ObjectCreationSystem::update(double deltaTime){
    std::cout << "update" << std::endl;
    m_ocean.fft_prime(m_time);
    Global::graphics.getShape("ocean")->updateVAO(m_ocean.get_vertices(), m_ocean.get_faces());
     m_time += m_timestep;

}
void ObjectCreationSystem::scrollEvent(double distance){}
void ObjectCreationSystem::mousePosEvent(double xpos, double ypos){}



