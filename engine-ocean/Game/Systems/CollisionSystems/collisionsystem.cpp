#include "collisionsystem.h"
#include "Game/Components/CollisionComponents/BoundingDynamicMesh.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include "Game/Components/TransformComponent.h"

CollisionSystem::CollisionSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects,
                                 std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects,
                                 std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects,
                                 std::map<int, Input>& input_map,
                                 std::map<std::string, BlackboardData>& global_blackboard) :
    m_gameobjects(gameobjects),
    m_dynamic_gameobjects(dynamic_gameobjects),
    m_rigid_gameobjects(rigid_gameobjects),
    m_uniform_grid_system(std::make_unique<UniformGrid>(dynamic_gameobjects)),
    m_global_blackboard(global_blackboard)
{

}

TransformComponent* CollisionSystem::getTransform(std::shared_ptr<GameObject> &go){
    return go->getComponent<TransformComponent>();

}

CollisionComponent* CollisionSystem::getCollisionComp(std::shared_ptr<GameObject> &go){
    return go->getComponent<CollisionComponent>();
}

glm::vec2 CollisionSystem::calculateCircleMVT(const Cylinder &a, const Cylinder &b){
    float len = glm::length(a.point - b.point);

    glm::vec2 mtv = ((b.point - a.point)/len) * (a.radius+b.radius - len);
    return mtv;
}

float CollisionSystem::calculateLineMVT(const Cylinder &a, const Cylinder &b){
    float aRight = b.max - a.min;
    float aLeft = a.max - b.min;
    if ((aLeft < 0) || (aRight < 0)){
        return -1.f;
    }
    if (aRight < aLeft){
        return aRight;
    }
    return -aLeft;
}

glm::vec3 CollisionSystem::collideCylinderCylinder(const std::shared_ptr<GameObject> &a_go, const std::shared_ptr<GameObject> &b_go){
    Cylinder a = a_go->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder();
    Cylinder b = b_go->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder();

    //std::cout << "a max: " << a.max << std::endl;
    //std::cout << "b max: " << b.max << std::endl;



    float len_squared = pow(a.point[0]-b.point[0],2) + pow(a.point[1]-b.point[1],2);
    float radius_sum = pow(a.radius + b.radius, 2);

    // check if circles overlap
    if (len_squared < radius_sum){
        //std::cout << "CIRCLES OVERLAP" << std::endl;
        // check if lines overlap
        if ((a.min < b.max) && (b.min < a.max)){
            std::cout << "COLLIDINGGG" << std::endl;
            // calculate both circle and line mtvs, and then return the shortest mtv
            glm::vec2 circle_mtv = calculateCircleMVT(a, b);
            float circle_len = glm::length(circle_mtv);

            float line_mtv = calculateLineMVT(a, b);

            if (abs(circle_len) < abs(line_mtv)){
                return glm::vec3(circle_mtv[0], 0.f, circle_mtv[1]);
            } else {
                return glm::vec3(0.f, line_mtv, 0.f);
            }
        }
    }

    // return translation of 0 if there is no collision
    return glm::vec3(0.f);
}


bool CollisionSystem::checkCollisionShape(const std::shared_ptr<GameObject> &go, std::string shape_type){
    // add more collision shapes here
    if (shape_type == "cylinder"){
        return (go->getComponent<CollisionComponent>()->hasCollisionShape<CylinderCollider>());
    }
    return false;
}

void CollisionSystem::resolveCollision(std::shared_ptr<GameObject> &go,
                                       glm::vec3 mtv,
                                       float deltaTime,
                                       std::string go_name){
    glm::vec3 potential_pos = go->getComponent<TransformComponent>()->getPos();
    go->getComponent<TransformComponent>()->old_pos = potential_pos;
    go->getComponent<TransformComponent>()->movingLaterally = false;

    // translate, and also update the collision cylinder point
    float a = go->getComponent<CollisionComponent>()->getAcceleration();
    float v = go->getComponent<CollisionComponent>()->getReboundVel();
    float delta_x = mtv.x + v*deltaTime + .5*a*pow(deltaTime, 2);
    float delta_z = mtv.z + v*deltaTime + .5*a*pow(deltaTime, 2);

    // update estimated_final_pos, so that tri-ellip collision can see if the new pos collides with environment
    glm::vec3 translationDir = glm::vec3(delta_x, mtv.y, delta_z);
    if (delta_x != 0 || delta_z != 0){
       go->getComponent<TransformComponent>()->movingLaterally = true;
    }

    std::shared_ptr<ModelTransform> temp_mt = go->getComponent<TransformComponent>()->getMT();
    temp_mt->translate(translationDir);
    potential_pos = temp_mt->getPos();


    m_global_blackboard[go_name].locationData.setToPos = potential_pos;
    go->getComponent<TransformComponent>()->estimated_final_pos = potential_pos;
}

void CollisionSystem::detectCylinderCollisions(std::shared_ptr<GameObject> &a,
                                               std::shared_ptr<GameObject> &b,
                                               float deltaTime,
                                               std::string a_name, std::string b_name){
    // if its a cyl-cyl collision
    glm::vec3 mtv(0.f);

    mtv = collideCylinderCylinder(a, b);
    if (mtv != glm::vec3(0.f)){
        resolveCollision(b, 1.f*mtv, deltaTime, b_name);
        resolveCollision(a, -1.f*mtv, deltaTime, a_name);
    }

}

bool CollisionSystem::areCollidable(const std::shared_ptr<GameObject> &a, const std::shared_ptr<GameObject> &b){
    return (a->hasComponent<CollisionComponent>() && b->hasComponent<CollisionComponent>());
}


// DYNAMIC-DYNAMIC COLLISIONS: CYLINDER
void CollisionSystem::dynamicDynamicCollisions(double deltaTime){
    std::set<std::set<std::string>> candidate_collison_set = m_uniform_grid_system->detectPossibleCollisions();

   // std::cout << "all collisions size : " << candidate_collison_set.size() << std::endl;

    for (const std::set<std::string> &set : candidate_collison_set){

       // std::cout << "collision group size : " << set.size() << std::endl;

        std::vector<std::string> coll_group;
        for (const std::string &name : set){
            coll_group.push_back(name);
        }

        for (int i=0; i < coll_group.size(); i++){
            auto a = m_dynamic_gameobjects.at(coll_group[i]);
               for (int j=i+1; j < coll_group.size(); j++){
                    auto b = m_dynamic_gameobjects.at(coll_group[j]);
                    // if a is not the same GO as b, and at least has collide components and thus is collidable,
                    if ((a != b) && (areCollidable(a, b))){
                        //std::cout << "collide: " << coll_group[i] << " + " << coll_group[j] << std::endl;
                        detectCylinderCollisions(a, b, deltaTime, coll_group[i], coll_group[j]);
                }
            }
        }
    }
}


void CollisionSystem::update(double deltaTime){
    // update uniform grid before detecting with it
    m_uniform_grid_system->updateAllGrid();
    dynamicDynamicCollisions(deltaTime);
}

void CollisionSystem::draw(){}
void CollisionSystem::scrollEvent(double distance){}
void CollisionSystem::mousePosEvent(double xpos, double ypos){}


