#include "bvhtree.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/boundingtriangle.h"
#include "Game/GameObjects/GameObject.h"
#include "Game/Systems/CollisionSystems/ellipsoidtrianglecollisionsystem.h"
#include <map>

BVHTree::BVHTree(std::map<std::string, std::shared_ptr<GameObject>>& rigid_gameobjects) :
    m_rigid_gameobjects(rigid_gameobjects)
{
    initializePrimitiveInfo();
    buildBVH();
    BVHNode one = m_bvhNodes[0];
    BVHNode two = m_bvhNodes[1];
    BVHNode three = m_bvhNodes[2];


   std::cout << "primitives size: " << m_primitives.size() << std::endl;
}

// populate primitive info vector -> primitives = triangles
void BVHTree::initializePrimitiveInfo(){
    for (const auto &go : m_rigid_gameobjects){
        N += go.second->getComponent<CollisionComponent>()->getCollisionShape<BoundingTriangle>()->getTriangleData().size();
        for (const Triangle &tri : go.second->getComponent<CollisionComponent>()->getCollisionShape<BoundingTriangle>()->getTriangleData()){
            // the order might not be accurate....
            BVHPrimitive prim = {tri};
            //std::cout << "centroid: (" << prim.centroid.x << ", " << prim.centroid.y << ", " << prim.centroid.z << ")" << std::endl;

            m_primitives.push_back(prim);

        }
    }
}

void BVHTree::buildBVH(){
    // initialize array of primitive indices
    m_primitive_indices = new int[N];
    for (int i=0; i<N; i++){
        m_primitive_indices[i] = i;
    }

    m_bvhNodes = new BVHNode[N*2 - 1];
    int rootNodeID = 0;
    m_nodesUsed = 1;

    // assign all primitives to root node
    BVHNode& root = m_bvhNodes[0];
    root.leftNode = 0;
    root.firstPrimitiveIndex = 0;
    root.primitiveCount = N;

    updateNodeBounds(rootNodeID);
    // subdivide recursively
    subdivide(rootNodeID);
}

void BVHTree::updateNodeBounds(int nodeIndex){
    BVHNode& node = m_bvhNodes[nodeIndex];
    node.bounds.max = glm::vec3(-INFINITY);
    node.bounds.min = glm::vec3(INFINITY);

    for (int first = node.firstPrimitiveIndex, i=0; i<node.primitiveCount; i++){
        int leafIndex = m_primitive_indices[first + i];
        BVHPrimitive& leafPrimitive = m_primitives[leafIndex];
        // update node bounds with the most max and most min coordinates
        node.bounds = getNodeBounds<Bounds3f>(node.bounds, leafPrimitive.bounds);
    }
}

void BVHTree::subdivide(int nodeIndex){
    BVHNode& node = m_bvhNodes[nodeIndex];
    // terminate recursion


    if (node.primitiveCount <= 1) return;

    // otherwise determine split axis and position
//    // SAHH HEREEEEE
//    int bestAxis = 0;
//    float bestPos = 0;
//    float bestCost = INFINITY;
//    for (int axis=0; axis<3; axis++){
//        for (int i=0; i<node.primitiveCount; i++){
//            BVHPrimitive& primitive = m_primitives[m_primitive_indices[node.firstPrimitiveIndex + i]];
//            float candidatePos = primitive.centroid[axis];
//            float cost = evaluateSAH(node, axis, candidatePos);
//            if (cost < bestCost){
//                bestPos = candidatePos;
//                bestAxis = axis;
//                bestCost = cost;
//            }

//        }
//    }

//    int axis = bestAxis;
//    float splitPos = bestPos;

//    glm::vec3 e = node.bounds.max - node.bounds.min;
//    float parentArea = (e.x*e.y) + (e.y*e.z) * (e.z*e.x);
//    float parentCost = node.primitiveCount*parentArea;
//    if (bestCost >= parentCost) return;

    // longest axis split (temporary)
    glm::vec3 extent = node.bounds.max - node.bounds.min;
    int axis = 0; // initialize to be x axis
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent[axis]) axis = 2;
    float splitPos = node.bounds.min[axis] + extent[axis] * 0.5f;

    // in-place partition
    int i = node.firstPrimitiveIndex;
    int j = i + node.primitiveCount - 1;
    while (i <= j){
        // if to the left of split axis
        if (m_primitives[m_primitive_indices[i]].centroid[axis] < splitPos){
            i++;
        } else {
            std::swap(m_primitive_indices[i], m_primitive_indices[j--]);
        }
    }

    // abort split if one of the sides is empty
    int leftCount = i - node.firstPrimitiveIndex;
    if (leftCount == 0 || leftCount == node.primitiveCount) return;

    // create child nodes
    int leftChildID = m_nodesUsed++;
    int rightChildID = m_nodesUsed++;

    // left child
    m_bvhNodes[leftChildID].firstPrimitiveIndex = node.firstPrimitiveIndex;
    m_bvhNodes[leftChildID].primitiveCount = leftCount;

    // right child
    m_bvhNodes[rightChildID].firstPrimitiveIndex = i;
    m_bvhNodes[rightChildID].primitiveCount = node.primitiveCount - leftCount;

    node.leftNode = leftChildID;
    node.primitiveCount = 0;

    // recurse
    updateNodeBounds(leftChildID);
    updateNodeBounds(rightChildID);
    subdivide(leftChildID);
    subdivide(rightChildID);
}

float BVHTree::evaluateSAH(BVHNode& node, int axis, float pos){
    BVHaabb leftBox, rightBox;
    leftBox.bounds.min = glm::vec3(INFINITY), rightBox.bounds.min = glm::vec3(INFINITY);
    leftBox.bounds.max = glm::vec3(-INFINITY), rightBox.bounds.max = glm::vec3(-INFINITY);

    int leftCount = 0, rightCount = 0;
    for (int i=0; i<node.primitiveCount; i++){
        BVHPrimitive& primitive = m_primitives[m_primitive_indices[node.firstPrimitiveIndex + i]];
        if (primitive.centroid[axis] < pos){
            leftCount++;
            leftBox.grow(getNodeBounds_Point<Bounds3f>(leftBox.bounds, primitive.triangle.vertexA));
            leftBox.grow(getNodeBounds_Point<Bounds3f>(leftBox.bounds, primitive.triangle.vertexB));
            leftBox.grow(getNodeBounds_Point<Bounds3f>(leftBox.bounds, primitive.triangle.vertexC));
        } else {
            rightCount++;
            rightBox.grow(getNodeBounds_Point<Bounds3f>(rightBox.bounds, primitive.triangle.vertexA));
            rightBox.grow(getNodeBounds_Point<Bounds3f>(rightBox.bounds, primitive.triangle.vertexB));
            rightBox.grow(getNodeBounds_Point<Bounds3f>(rightBox.bounds, primitive.triangle.vertexC));
        }
    }
    float cost = leftCount * leftBox.area() + rightCount * rightBox.area();
    return cost > 0 ? cost : INFINITY;
}


void BVHTree::intersectBVH(const glm::vec3 posA, const glm::vec3 posB, const int nodeIndex, std::vector<Triangle> &candidates, glm::vec3 ellip_R){
    //std::cout << "-" << std::endl;

    BVHNode& node = m_bvhNodes[nodeIndex];
    if (!intersectsNode(posA, posB, node.bounds.min, node.bounds.max, ellip_R+glm::vec3(.001))){
        return;
    }
    if (node.isLeaf()){
        // intersect with primitives of the leaf
        //std::cout << "candidate size: " << node.primitiveCount << std::endl;
        for (int i=0; i<node.primitiveCount; i++){
            Triangle triangle = m_primitives[m_primitive_indices[node.firstPrimitiveIndex + i]].triangle;
//            BVHPrimitive prim = m_primitives[m_primitive_indices[node.firstPrimitiveIndex + i]];

//            if (prim.centroid == glm::vec3(39.42457, 1.84954, -15.7088)){
//                std::cout << "found correct triangle" << std::endl;
//            }
            //std::cout << "centroid: (" << prim.centroid.x << ", " << prim.centroid.y << ", " << prim.centroid.z << ")" << std::endl;
            candidates.push_back(triangle);

        }
    } else {
        // recurse
        intersectBVH(posA, posB, node.leftNode, candidates, ellip_R);
        intersectBVH(posA, posB, node.leftNode + 1, candidates, ellip_R);
    }
}

std::vector<Triangle> BVHTree::getBVHDetectedCollisions(glm::vec3 posA, glm::vec3 posB, glm::vec3 ellip_R){
    // keep recursing, starting at root node, and then return final result in candidate_obstacles
    std::vector<Triangle> candidate_rigid_tri;
    intersectBVH(posA, posB, 0, candidate_rigid_tri, ellip_R);
    return candidate_rigid_tri;
}

// determines if movement ray intersects an AABB via the slab test
bool BVHTree::intersectsNode(glm::vec3 posA, glm::vec3 posB, glm::vec3 min, glm::vec3 max, glm::vec3 ellip_R){
    glm::vec3 object_min = glm::vec3(std::min(posA.x, posB.x), std::min(posA.y, posB.y), std::min(posA.z, posB.z));
    glm::vec3 object_max = glm::vec3(std::max(posA.x, posB.x), std::max(posA.y, posB.y), std::max(posA.z, posB.z));
    object_min -= ellip_R;
    object_max += ellip_R;

    bool x_int = object_max.x > min.x && object_min.x < max.x;
    bool y_int = object_max.y > min.y && object_min.y < max.y;
    bool z_int = object_max.z > min.z && object_min.z < max.z;

    if (x_int && y_int && z_int){
        return true;
    }
    return false;
}

BVHTree::~BVHTree(){
    delete[] m_bvhNodes;
    delete[] m_primitive_indices;
    m_bvhNodes = NULL;
    m_primitive_indices = NULL;

}

