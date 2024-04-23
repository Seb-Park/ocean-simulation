#include "boundingdynamicmesh.h"
#include "Graphics/modeltransform.h"
#include <vector>

// without mesh obj data
BoundingDynamicMesh::BoundingDynamicMesh(std::shared_ptr<ModelTransform> mt,
                                         const glm::vec3 &initial_pos) :
    m_mt(mt)
{
    m_isMesh = false;
    m_bounding_cylinder = std::make_shared<CylinderCollider>(initial_pos, mt->getScale());
}

// with mesh obj data
BoundingDynamicMesh::BoundingDynamicMesh(std::shared_ptr<ModelTransform> mt,
                                         const glm::vec3 &initial_pos,
                                         std::vector<glm::vec3> &obj_data) :
    m_mt(mt),
    m_obj_data(obj_data)
{
    m_isMesh = true;
    m_bounding_cylinder = std::make_shared<CylinderCollider>(initial_pos, mt->getScale()*getMeshDimensions()*2.f);
}


glm::vec3 BoundingDynamicMesh::getCenterPos(){
    return m_mt->getPos();
}

void BoundingDynamicMesh::updateCenterPos(glm::vec3 new_pos){
    m_mt->setPos(new_pos);
    m_bounding_cylinder->updateCollisionPos(new_pos);
}



glm::vec3 BoundingDynamicMesh::getMeshDimensions(){
    float max_x =  m_obj_data[0].x;
    float min_x = m_obj_data[0].x;
    float max_y = m_obj_data[0].y;
    float min_y = m_obj_data[0].y;
    float max_z = m_obj_data[0].z;
    float min_z = m_obj_data[0].z;

    for (const glm::vec3 &v : m_obj_data){
        // check max
        if (v.x > max_x){
            max_x = v.x;
        }
        if (v.y > max_y){
            max_y = v.y;
        }
        if (v.z > max_z){
            max_z = v.z;
        }

        // check mins
        if (v.x < min_x){
            min_x = v.x;
        }
        if (v.y < min_y){
            min_y = v.y;
        }
        if (v.z < min_z){
            min_z = v.z;
        }
    }

    float r_x = (max_x - min_x)/2.f;
    float r_y = (max_y - min_y)/2.f;
    float r_z = (max_z - min_z)/2.f;

    return glm::vec3(r_x, r_y, r_z);
}

glm::vec3 BoundingDynamicMesh::getEllipsoidDimensions(){
    if (m_isMesh){
        return m_mt->getScale()*getMeshDimensions();
    }

    return m_mt->getScale()/2.f;
}

Cylinder BoundingDynamicMesh::getCylinder(){
    return m_bounding_cylinder->getCylinder();
}


