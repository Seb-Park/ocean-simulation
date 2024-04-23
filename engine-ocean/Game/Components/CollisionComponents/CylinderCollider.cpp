#include "cylindercollider.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Graphics/modeltransform.h"

CylinderCollider::CylinderCollider(glm::vec3 initial_pos, glm::vec3 initial_scale){
    m_scale = initial_scale;

    m_cyl.point = glm::vec2(initial_pos.x, initial_pos.z);
    m_cyl.radius = .5f * glm::max(initial_scale.x, initial_scale.z);
    m_cyl.height = 1.f * abs(initial_scale.y);
    m_cyl.min = initial_pos.y - (initial_scale.y/2.f);
    m_cyl.max = initial_pos.y + (initial_scale.y/2.f);

    m_cyl.aabbDimensions = abs(glm::vec3(2*m_cyl.radius, m_cyl.height, 2*m_cyl.radius));
    m_cyl.aabbCenterPos = initial_pos;

    updateBounds();
}

void CylinderCollider::updateCollisionPos(glm::vec3 new_pos){
    m_cyl.point = glm::vec2(new_pos.x, new_pos.z); // x and z loc
    m_cyl.radius = .5f * m_scale.x; // x and z dim
    m_cyl.height = 1.f * m_scale.y; // y dimensions
    m_cyl.min = new_pos.y - (m_scale.y/2.f); // y coord
    m_cyl.max = new_pos.y + (m_scale.y/2.f); // y coord
    m_cyl.aabbCenterPos = new_pos;

    updateBounds();
}


void CylinderCollider::updateBounds(){
    Bounds3f bounds;
    bounds.min = glm::vec3(m_cyl.aabbCenterPos-glm::vec3(m_cyl.aabbDimensions/2.f));
    bounds.max = glm::vec3(m_cyl.aabbCenterPos+glm::vec3(m_cyl.aabbDimensions/2.f));

    m_cyl.bounds = bounds;
}


Cylinder CylinderCollider::getCylinder(){
    return m_cyl;
}




