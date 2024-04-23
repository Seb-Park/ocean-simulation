#include "boundingtriangle.h"
#include "Graphics/global.h"
#include "Graphics/modeltransform.h"
#include <iostream>
#include "glm/gtx/hash.hpp"

// ONLY FOR ENVIRONMENTS
BoundingTriangle::BoundingTriangle(const std::vector<glm::vec3> &obj_data,
                                   const std::shared_ptr<ModelTransform> &mt,
                                   bool isGround) :
    obj_mt(mt),
    m_isGround(isGround)
{
    populateTriangleData(obj_data);
    calculateBounds(obj_data);
    //m_datasize = obj_data.size();
}

glm::vec3 BoundingTriangle::getRandomSurfacePos(){
    int randomIndex = std::floor(Global::graphics.generateRandomNumbers(0, m_triangles.size()-1));
    int randomVertex = std::floor(Global::graphics.generateRandomNumbers(0, 3));

    Triangle randomTri = m_triangles[randomIndex];

    if (glm::dot(glm::vec3(0,1,0), randomTri.normal) > 0.2f){
        switch(randomVertex){
            case 0:
                return randomTri.vertexA;
            break;
            case 1:
                return randomTri.vertexB;
            break;
            default:
                return randomTri.vertexC;
            break;
        }
    } else {
        // do again until returning a surface triangle
        getRandomSurfacePos();
    }
}


void BoundingTriangle::addTriangle(const glm::vec3 &vertexA, const glm::vec3 &vertexB, const glm::vec3 &vertexC){
    Triangle tri;
    tri.vertexA = vertexA;
    tri.vertexB = vertexB;
    tri.vertexC = vertexC;

    tri.edge1 = vertexB - vertexA; // edge ab
    tri.edge2 = vertexC - vertexA; // edge bc
    tri.normal = glm::normalize(glm::cross(tri.edge1, tri.edge2));
    tri.bounds = calculateTriangleBounds(vertexA, vertexB, vertexC);

    // if triangle is a ground triangle
    if (m_isGround){
       if (glm::dot(glm::vec3(0,1,0), tri.normal) > 0.2f){
           //std::cout << "area: " << getArea(vertexA, vertexB, vertexC) << std::endl;
           tesselateTriangle(vertexA, vertexB, vertexC);
       }
    }

    m_triangles.push_back(tri);
}

float BoundingTriangle::getArea(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C){
    glm::vec3 AB = B-A;
    glm::vec3 AC = C-A;
    return .5f*glm::length((glm::cross(AB,AC)));

}

glm::vec3 BoundingTriangle::getCentroid(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C){
    return .333f*(A + B + C);
}

// tesselation into smaller triangles
void BoundingTriangle::tesselateTriangle(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C){
    float min_area = 1.5f;

    // add centroid if area is small enough
    if (getArea(A,B,C) <= min_area){
        m_surface_points.push_back(getCentroid(A, B, C));
        return;
    }

    // otherwise divide triangle in 4 (tesselate)

        glm::vec3 ab_mid = .5f*(A+B);
        glm::vec3 ac_mid = .5f*(A+C);
        glm::vec3 bc_mid = .5f*(B+C);
        m_surface_points.push_back(ab_mid);
        m_surface_points.push_back(ac_mid);
        m_surface_points.push_back(bc_mid);

        tesselateTriangle(ab_mid, ac_mid, A);
        tesselateTriangle(bc_mid, ab_mid, B);
        tesselateTriangle(bc_mid, ac_mid, C);
        tesselateTriangle(ab_mid, ac_mid, bc_mid);
}


void BoundingTriangle::populateTriangleData(const std::vector<glm::vec3> &obj_data){
    for (int i=0; i<obj_data.size(); i += 3){
        // convert to worldspace
        glm::mat4 modelMat = obj_mt->getModelMatrix();
        glm::vec3 v1 = modelMat*glm::vec4(obj_data[i],1.0);
        glm::vec3 v2 = modelMat*glm::vec4(obj_data[i+1],1.0);
        glm::vec3 v3 = modelMat*glm::vec4(obj_data[i+2],1.0);

        // make triangle
        addTriangle(v1, v2, v3);
    }

//    m_surface_points.reserve(m_unique_surface_points.size());
//    std::copy(m_unique_surface_points.begin(), m_unique_surface_points.end(), m_surface_points.begin());
}

Bounds3f BoundingTriangle::calculateTriangleBounds(const glm::vec3 &vertexA,
                                                   const glm::vec3 &vertexB,
                                                   const glm::vec3 &vertexC){
    Bounds3f bounds;

    float max_x, min_x = vertexA.x;
    float max_y, min_y = vertexA.y;
    float max_z, min_z = vertexA.z;

    if (vertexB.x > max_x) max_x = vertexB.x;
    if (vertexB.y > max_y) max_y = vertexB.y;
    if (vertexB.z > max_z) max_z = vertexB.z;

    if (vertexC.x > max_x) max_x = vertexC.x;
    if (vertexC.y > max_y) max_y = vertexC.y;
    if (vertexC.z > max_z) max_z = vertexC.z;

    if (vertexB.x < min_x) min_x = vertexB.x;
    if (vertexB.y < min_y) min_y = vertexB.y;
    if (vertexB.z < min_z) min_z = vertexB.z;

    if (vertexC.x < min_x) min_x = vertexC.x;
    if (vertexC.y < min_y) min_y = vertexC.y;
    if (vertexC.z < min_z) min_z = vertexC.z;


    bounds.max = glm::vec3(max_x, max_y, max_z);
    bounds.min = glm::vec3(min_x, min_y, min_z);

    return bounds;
}

void BoundingTriangle::calculateBounds(const std::vector<glm::vec3> &obj_data){
    max_x = obj_data[0].x;
    min_x = obj_data[0].x;
    max_y = obj_data[0].y;
    min_y = obj_data[0].y;
    max_z = obj_data[0].z;
    min_z = obj_data[0].z;

    for (const glm::vec3 &v : obj_data){
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
}

Bounds3f BoundingTriangle::getMeshBounds(){
    Bounds3f bounds;
    bounds.max = glm::vec3(max_x, max_y, max_z);
    bounds.min = glm::vec3(min_x, min_y, min_z);
    return bounds;
}

std::vector<glm::vec3> BoundingTriangle::getSurfacePoints(){
    if (!m_isGround){
        std::cout << "getting surface points of not-ground object!" << std::endl;
    }
    return m_surface_points;
}

std::vector<Triangle> BoundingTriangle::getTriangleData(){
    return m_triangles;
}
