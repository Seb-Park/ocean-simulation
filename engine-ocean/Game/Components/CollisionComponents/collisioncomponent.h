#ifndef COLLISIONCOMPONENT_H
#define COLLISIONCOMPONENT_H

#include "Game/Components/CollisionComponents/BoundingShape.h"
#include <map>
#include <string>
#include "Game/TypeMap.h"
#include "Graphics/modeltransform.h"
#include "glm/glm.hpp"
#include "Game/Components/Component.h"


class CollisionComponent : public Component
{
public:
    template <typename T>
    void addCollisionShape(std::shared_ptr<T> &&component){
        m_collision_shapes.put<T>(std::forward<std::shared_ptr<T>>(component));
    }

    template <typename T>
    bool hasCollisionShape(){
        return m_collision_shapes.contains<T>();
    }

    template <class T>
    T* getCollisionShape(){
        auto comp = m_collision_shapes.find<T>();
        assert(comp != m_collision_shapes.end());
        return static_cast<T*>(comp->second.get());
    }

    CollisionComponent(std::string shapeType, std::shared_ptr<ModelTransform> mt, const glm::vec3 &initial_pos, std::vector<glm::vec3> &obj_data);
    CollisionComponent(std::string shapeType, std::shared_ptr<ModelTransform> mt, const glm::vec3 &initial_pos);
    CollisionComponent(std::string shapeType, const std::vector<glm::vec3> &obj_data, const std::shared_ptr<ModelTransform> &mt,
                       bool isGround = false);

    bool isRigidBody();
    float getReboundVel();
    float getAcceleration();

private:
    TypeMap<std::shared_ptr<BoundingShape>> m_collision_shapes;
    bool m_isRigid = true;
    float m_rebound_vel = 0.f;
    float m_acceleration = 0.f;
};

#endif // COLLISIONCOMPONENT_H
