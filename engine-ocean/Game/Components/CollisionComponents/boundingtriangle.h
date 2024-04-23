#ifndef BOUNDINGTRIANGLE_H
#define BOUNDINGTRIANGLE_H
#include "Graphics/modeltransform.h"
#include "glm/glm.hpp"
#include <set>
#include <vector>
#include "BoundingShape.h"

struct Bounds3f {
    glm::vec3 min;
    glm::vec3 max;
};

struct Triangle{
    glm::vec3 vertexA; // one vertex
    glm::vec3 vertexB;
    glm::vec3 vertexC;

    glm::vec3 edge1; // two edges
    glm::vec3 edge2;

    glm::vec3 normal;// = glm::cross(edge1, edge2);
    Bounds3f bounds;


};

class BoundingTriangle : public BoundingShape
{
public:
    BoundingTriangle(const std::vector<glm::vec3> &obj_data,
                     const std::shared_ptr<ModelTransform> &mt,
                     bool isGround = false);
    std::vector<Triangle> getTriangleData();

    Bounds3f getMeshBounds();
    std::vector<glm::vec3> getSurfacePoints();
    glm::vec3 getRandomSurfacePos();



private:
    void addTriangle(const glm::vec3 &vertexA, const glm::vec3 &vertexB, const glm::vec3 &vertexC);
    void populateTriangleData(const std::vector<glm::vec3> &obj_data);
    void calculateBounds(const std::vector<glm::vec3> &obj_data);
    Bounds3f calculateTriangleBounds(const glm::vec3 &vertexA, const glm::vec3 &vertexB, const glm::vec3 &vertexC);

    float getArea(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C);
    glm::vec3 getCentroid(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C);
    void tesselateTriangle(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C);


    int m_datasize = 0;





    std::vector<Triangle> m_triangles;
    std::shared_ptr<ModelTransform> obj_mt;
    std::set<glm::vec3> m_unique_surface_points;
    std::vector<glm::vec3> m_surface_points;
    bool m_isGround = false;

    float min_x, min_y, min_z, max_x, max_y, max_z;
};

#endif // BOUNDINGTRIANGLE_H
