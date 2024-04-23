#include "collisioncomponent.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include "Game/Components/CollisionComponents/BoundingDynamicMesh.h"
#include "Graphics/shape.h"
#include <memory>

// for dynamic objects
CollisionComponent::CollisionComponent(std::string shapeType, std::shared_ptr<ModelTransform> mt, const glm::vec3 &initial_pos)
{
    if (shapeType == "dynamic_mesh"){
        //addCollisionShape<CylinderCollider>(std::make_shared<CylinderCollider>(initial_pos, initial_scale));
        addCollisionShape<BoundingDynamicMesh>(std::make_shared<BoundingDynamicMesh>(mt, initial_pos));
    }
}

CollisionComponent::CollisionComponent(std::string shapeType, std::shared_ptr<ModelTransform> mt, const glm::vec3 &initial_pos, std::vector<glm::vec3> &obj_data)
{
    if (shapeType == "dynamic_mesh"){
        addCollisionShape<BoundingDynamicMesh>(std::make_shared<BoundingDynamicMesh>(mt, initial_pos, obj_data));
    }
}

// for rigid meshes / environment
CollisionComponent::CollisionComponent(std::string shapeType,
                                       const std::vector<glm::vec3> &obj_data,
                                       const std::shared_ptr<ModelTransform> &mt,
                                       bool isGround){
    if (shapeType == "obj"){
        addCollisionShape<BoundingTriangle>(std::make_shared<BoundingTriangle>(obj_data, mt, isGround));
    }

}



bool CollisionComponent::isRigidBody(){
    return m_isRigid;
}

float CollisionComponent::getReboundVel(){
    return m_rebound_vel;
}

float CollisionComponent::getAcceleration(){
    return m_acceleration;
}
