#ifndef ELLIPSOIDTRIANGLECOLLISIONSYSTEM_H
#define ELLIPSOIDTRIANGLECOLLISIONSYSTEM_H
#include "Game/Components/CollisionComponents/BoundingEllipsoid.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/Systems/CollisionSystems/BVH/bvhtree.h"
#include "glm/glm.hpp"
#include <vector>

struct CollisionData{
    bool hasCollided = true; // only false when there is no collision data to return
    bool isVertex = false;
    float t = INFINITY;
    glm::vec3 triangle_n; // normal of the triangle its colliding with
    glm::vec3 intersection_pos;
    glm::vec3 point_of_contact;
};

class EllipsoidTriangleCollisionSystem
{
public:
    EllipsoidTriangleCollisionSystem(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects);


    std::pair<std::vector<CollisionData>, glm::vec3> mtvSlide(glm::vec3 initial_pos, glm::vec3 final_pos, glm::vec3 dimensions);

private:
    CollisionData collideWithWorld(const glm::vec3 &posA, const glm::vec3 &posB);
    void setEllipsoidDimensions(const glm::vec3 &dimensions);
    glm::vec3 doNudge(glm::vec3 curr_pos, CollisionData collision);
    CollisionData revertDataToOriginalSpace(const CollisionData &sphere_datum, const glm::vec3 &ellip_R, const glm::vec3 orig_triangle_n);
    CollisionData detectTriangleCollision(const Triangle &sphere_tri, const glm::vec3 &posA, const glm::vec3 &posB);
    CollisionData getNearestCollision(const std::vector<CollisionData> &collisionData, const Triangle &sphere_tri);
    glm::mat3 makeConversionMT(glm::vec3 initial_center_pos, glm::vec3 dimensions);


    bool detectSphere_with_Vertex(const glm::vec3 &v, const glm::vec3 triangle_normal, const glm::vec3 &posA, const glm::vec3 &posB);
    bool detectSphere_with_Edge(const glm::vec3 &v1, const glm::vec3 &v2,
                                const glm::vec3 &triangle_normal,
                                 const glm::vec3 &posA, const glm::vec3 &posB);
    void makeCollisionData(const float t, const glm::vec3 triangle_n, const glm::vec3 intersection_pos, const glm::vec3 point_of_contact, const bool isVertex);
    bool detectSphere_with_Interior(const Triangle &triangle, const glm::vec3 posA, const glm::vec3 posB);


    bool isSameDirection(const glm::vec3 edge1, const glm::vec3 edge2, const glm::vec3 n);
    glm::vec3 getTriangleEdge(const glm::vec3 v1, const glm::vec3 v2);
    float raycastSphere(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const glm::vec3 posA, const glm::vec3 posB);
    float raycastTriangle(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const Triangle &triangle);
    float solveQuadratic(float a, float b, float c);
    glm::vec3 getIntersectionPt(const glm::vec3 &ray_p, const glm::vec3 &ray_d, const float t);
    Triangle convertTriangleToSphereSpace(const glm::vec3 &orig_R, const Triangle &orig_triangle);

    std::vector<CollisionData> m_collisionData;
    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;

    glm::vec3 m_ellipsoid_dimensions = glm::vec3(1.0);

    float ZERO = -.0001f;
    float ONE = 1.0001f;

    std::shared_ptr<ModelTransform> m_conversion_mt;
    glm::mat3 m_scale_normals_mt = glm::mat3(1.f);

    std::vector<std::shared_ptr<GameObject>> m_candidate_obstacles;
    std::unique_ptr<BVHTree> m_bvh_system;

};

#endif // ELLIPSOIDTRIANGLECOLLISIONSYSTEM_H
