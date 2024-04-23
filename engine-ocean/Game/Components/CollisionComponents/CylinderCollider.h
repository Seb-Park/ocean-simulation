#ifndef CYLINDERCOLLIDER_H
#define CYLINDERCOLLIDER_H

#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Graphics/modeltransform.h"
#include "glm/glm.hpp"
#include "BoundingShape.h"
#include <memory>


struct Cylinder {
    glm::vec2 point; // bottom Center
    float radius;

    // lines
    float height;
    float min;
    float max;

    glm::vec3 aabbDimensions;
    glm::vec3 aabbCenterPos;
    Bounds3f bounds;
};

class CylinderCollider
{
public:
    CylinderCollider(glm::vec3 initial_pos, glm::vec3 initial_scale);
    Cylinder getCylinder();// override;
    void updateCollisionPos(glm::vec3 new_pos);
    void updateBounds();
private:
    glm::vec3 m_scale;
    Cylinder m_cyl;
};

#endif // CYLINDERCOLLIDER_H
