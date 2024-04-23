#ifndef BVHTREE_H
#define BVHTREE_H
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/GameObjects/GameObject.h"
#include "glm/glm.hpp"
#include <map>
#include <vector>

// primitive = triangle
struct BVHPrimitive{

    Bounds3f bounds;
    glm::vec3 centroid;
    Triangle triangle;

    BVHPrimitive(const Triangle &tri):
        triangle(tri),
        bounds(tri.bounds),
        centroid(.333f*(tri.vertexA + tri.vertexB + tri.vertexC)){}
};

struct BVHNode{
   // glm::vec3 aabbMin, aabbMax;
    Bounds3f bounds;
    int leftNode, firstPrimitiveIndex, primitiveCount;
    bool isLeaf(){
        return primitiveCount > 0;
    }
};

struct BVHaabb {
    Bounds3f bounds;


    void grow(Bounds3f newbounds){
        bounds.min = newbounds.min;
        bounds.max = newbounds.max;
    }

    float area(){
        glm::vec3 e = bounds.max-bounds.min;
        return (e.x*e.y) + (e.y*e.z) * (e.z*e.x);
    }
};


template <typename T> Bounds3f
getNodeBounds(const Bounds3f &b1, const Bounds3f &b2){
    Bounds3f bounds;
    bounds.min = glm::vec3(std::min(b1.min.x, b2.min.x),
                           std::min(b1.min.y, b2.min.y),
                           std::min(b1.min.z, b2.min.z));
    bounds.max = glm::vec3(std::max(b1.max.x, b2.max.x),
                           std::max(b1.max.y, b2.max.y),
                           std::max(b1.max.z, b2.max.z));
    return bounds;
}

template <typename T> Bounds3f
getNodeBounds_Point(const Bounds3f &b, const glm::vec3 &p){
    Bounds3f bounds;
    bounds.min = glm::vec3(std::min(b.min.x, p.x),
                           std::min(b.min.y, p.y),
                           std::min(b.min.z, p.z));
    bounds.max = glm::vec3(std::max(b.max.x, p.x),
                           std::max(b.max.y, p.y),
                           std::max(b.max.z, p.z));
    return bounds;
}

class BVHTree

{
    BVHNode *m_bvhNodes;
    int *m_primitive_indices;
public:
    BVHTree(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects);
    ~BVHTree();
    std::vector<Triangle> getBVHDetectedCollisions(glm::vec3 posA, glm::vec3 posB, glm::vec3 ellip_R);
    void intersectBVH(glm::vec3 posA, glm::vec3 posB, const int nodeIndex, std::vector<Triangle> &candidates, glm::vec3 ellip_R);


private:
    float evaluateSAH(BVHNode& node, int axis, float pos);
    void initializePrimitiveInfo();
    void buildBVH();
    void updateNodeBounds(int nodeIndex);
    void subdivide(int nodeIndex);
    bool intersectsNode(glm::vec3 posA, glm::vec3 posB, glm::vec3 min, glm::vec3 max, glm::vec3 ellip_R);


    //std::unique_ptr<EllipsoidTriangleCollisionSystem> m_ellipsoid_triangle_collision_system;

    std::vector<BVHPrimitive> m_primitives;
    std::map<std::string, std::shared_ptr<GameObject>>& m_rigid_gameobjects;


    //std::vector<BVHNode> m_bvhNodes;
    //std::vector<int> m_primitive_indices;
    int m_nodesUsed = 1;
    int N = 0;

   // std::pair<std::vector<CollisionData>, glm::vec3> m_detected_collisions;
};

#endif // BVHTREE_H
