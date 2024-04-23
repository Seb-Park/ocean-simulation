#ifndef BOUNDINGDYNAMICMESH_H
#define BOUNDINGDYNAMICMESH_H
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include "Graphics/modeltransform.h"
#include "glm/glm.hpp"
#include <memory>
#include <vector>


class BoundingDynamicMesh : public BoundingShape
{
public:
    BoundingDynamicMesh(std::shared_ptr<ModelTransform> mt,
                                             const glm::vec3 &initial_pos);
    BoundingDynamicMesh(std::shared_ptr<ModelTransform> mt,
                                             const glm::vec3 &initial_pos,
                                             std::vector<glm::vec3> &obj_data);
    glm::vec3 getCenterPos();
    void updateCenterPos(glm::vec3 new_pos);
    glm::vec3 getEllipsoidDimensions();
    Cylinder getCylinder();

private:
    glm::vec3 getMeshDimensions();

    bool m_isMesh = false;
    std::shared_ptr<ModelTransform> m_mt;
    std::vector<glm::vec3> m_obj_data;
    std::shared_ptr<CylinderCollider> m_bounding_cylinder;
};

#endif // BOUNDINGDYNAMICMESH_H
